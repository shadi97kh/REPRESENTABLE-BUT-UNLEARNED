# Attainment from the five observed groups

**PROVED HERE:** for every fixed known \(\theta\in[1/2,2]\) and \(0<\delta\le1/3{,}200{,}000\), and every fixed truth in the full original profile and coefficient class in [problem.md](problem.md), the statistic below satisfies
\[
\widehat I_N-I_\eta
=\frac1N\sum_{a=0}^4\sum_{i=1}^n\psi_{\eta,a}(Y_{ai})
+o_{P_\eta}(N^{-1/2}),\qquad N=5n.
\tag{1}
\]
Each source component has finite second moment. The expansion is regular along every fixed admissible differentiable-in-quadratic-mean path specified below. Neither a profile neighborhood shrinking with separation nor an oracle coefficient root is used. The quantitative norm bound is [conditioning.md](conditioning.md), (10).

This is a pointwise asymptotic and mathematical arithmetic construction. It is not an MSE, efficiency, uniform shrinking-separation, fixed-resolution-observation or practical-runtime theorem. No estimator, numerical integration or sampled-data calculation was executed for this document.

## 1. Observable specification and simultaneous schedules

Write \(J=[-2,0]\), \(r=-1\), and distinguish the fixed calibration radius \(R\) from the growing target cutoff \(B_N\). All definitions of \(C,E,D_K,H_m,G_{1,K},G_{2,K}\), the calibration map and its fixed transport charts are in [commutator_and_calibration.md](commutator_and_calibration.md). Use those exact signs.

The population constants are public:
\[
\lambda=1-\kappa\delta/2,\quad \kappa=1/64,\quad
\gamma=\tfrac12e^{.4}-2e^{-1.9},\quad
\tau=(2e^{1.1})^{-1},\quad q=1-\tau\gamma.
\]
Take
\[
K=\left\lceil\frac{2\log N}{|\log\lambda|}\right\rceil+2,\quad
B_N=\sqrt{1.5\log N},\quad M_N=\lceil4B_N+30\rceil+2,\quad
k_N=\left\lceil\frac{\log N}{|\log q|}\right\rceil .
\tag{2}
\]
In particular \(K\le\lceil4\log N/(\kappa\delta)\rceil+2\). The notation \(K=O(\log N)\) henceforth always fixes \(\delta,\theta\).

For each group sort the observed positive responses and define \(\widetilde Q_a\) by linear interpolation through \((i/n,Y_{a(i)})\), constant on \([0,1/n]\). Put \(\widetilde F_a(z)=\widetilde Q_a(\Phi(z))\). The latent values used in probability-integral-transform arguments below are never inputs.

Choose the four calibration transport integers once from the public query intervals, with the rational-enclosure convention in the companion proof. More explicitly, enclose each left endpoint in an interval of width at most \(1/100\), apply the integer ceiling to its rational lower endpoint, and keep this integer for the whole coefficient box. The resulting translated intervals lie in \([-7/4,-1/25]\subset J\). Use these same charts in the mathematical statistic and its numerical realization. Exact comparison of an arbitrary computable real to an integer is not required.

For each pair start \(c_i^0=(.75,.75)\). Form
\[
\widehat y_i=(e^R\widetilde F_i(-R),e^{-h(R)}\widetilde F_i(R))^\top
\]
and form \(\widehat B_i(c)\) by replacing \(g_1,g_2\) in the scaled map \(B\) with \(G_{1,K}(\widetilde F),G_{2,K}(\widetilde F)\) using the fixed charts. Apply exactly \(2k_N\) public updates
\[
c_i^{t+1}=\Pi_{[.4,1.1]^2}
\{c_i^t+\tau[\widehat y_i-\widehat B_i(c_i^t)]\},\qquad
\widehat c_i=c_i^{2k_N}.
\tag{3}
\]
This is a finite data-defined recursion. It uses no true profile derivative, inverse true Jacobian, fitted profile, nearest-oracle-root rule or target observations.

Return zero for \(n<2\), invalid nonpositive data or the observable cap failure \(\max Y>N^2\). Otherwise compute the finite target (7) and clip it to \([-L_I,L_I]\), where the public constant
\[
L_I=1+2(e+1)^2(m_0+\overline m_{h,1}),\qquad m_0=e^{1/2},
\tag{4}
\]
uses the moment bound in conditioning (2). This clipping level exceeds the true target magnitude throughout the original class. The finite-bit gate is specified separately in Section 8.

## 2. Exact signed target and finite remainder

**PROVED HERE.** Let \(\rho_1=\phi\), \(\rho_2(x)=\phi(h^{-1}x)/h'(h^{-1}x)\), and
\[
K_{1,c}(x)=\rho_1(x-\alpha_1-\alpha_2)-\rho_1(x-\alpha_1)
-\rho_1(x-\alpha_2)+\rho_1(x),
\]
with \(K_{2,c}\) defined in the same way using \(\rho_2,\beta_1,\beta_2\).
Changing variables in the four expectations gives
\[
I(c,g)=\int K_{1,c}g_1+\int K_{2,c}g_2
=\int A_c(z)F_0(z)\,dz+\int\mathcal L_c(z)g_1(z)\,dz,
\]
where \(A_c=h'K_{2,c}\circ h\), \(\mathcal L_c=K_{1,c}-A_c\).
All integrals are absolutely convergent by Section 5. Both signed kernels, and hence \(\mathcal L_c\), have integral zero.

Choose the unit cell \([-1,0]\subset J\), and set
\[
P_c(t)=\sum_{j\in\mathbb Z}\mathcal L_c(t+j),\qquad
W_c(z)=
\begin{cases}
\sum_{j\ge1}\mathcal L_c(z+j),&z\ge-1,\\
-\sum_{j\ge0}\mathcal L_c(z-j),&z<-1.
\end{cases}
\tag{5}
\]
Partition the line into translates of this cell, and telescope
\(g_1(x+1)-g_1(x)=F_3(x)-F_0(x)\) outward. Absolute convergence permits reordering. Since \(\int_{-1}^0P_c=0\), the result is
\[
I=\int_{\mathbb R}\{(A_c-W_c)F_0+W_cF_3\}\,dz
+\int_{-1}^0P_c(t)[g_1(t)-g_1(r)]\,dt.
\tag{6}
\]
The join is \(W_c(-1+)-W_c(-1-)=P_c(-1)\), not zero in general.

Let \(P_{M,c}\) sum over \(-M,\ldots,M\), and \(W_{M,c}\) over \(1,\ldots,M\) on the positive branch and \(0,\ldots,M-1\) on the negative branch. Define
\[
H_N(\widetilde F,c)=
\int_{-B_N}^{B_N}\{(A_c-W_{M,c})\widetilde F_0+
W_{M,c}\widetilde F_3\}\,dz
+\int_{-1}^0P_{M,c}(t)D_K(\widetilde F;t,r)\,dt.
\tag{7}
\]
Its integrals are initially exact. The estimator uses \(c=\widehat c\).
The finite-periodization constant is retained: generally \(\int P_{M,c}\ne0\), so the central integrand is the profile difference \(D_K\), never an unnormalized profile with its constant dropped.

For clarity, the exact difference \(I-H_N(F,c)\) is the sum of
\[
\begin{split}
&\int_{|z|>B_N}\{(A_c-W_c)F_0+W_cF_3\}\,dz,\\
&\int_{-B_N}^{B_N}(W_c-W_{M,c})(F_3-F_0)\,dz,\\
&\int_{-1}^0(P_c-P_{M,c})[g_1(t)-g_1(r)]\,dt,\\
&\int_{-1}^0P_{M,c}(t)[g_1(C^Kt)-g_1(C^Kr)]\,dt.
\end{split}
\tag{8}
\]
Thus signs, terminal terms, the shared reference and the join survive finite truncation explicitly.

## 3. Compact empirical quantiles and growing series

**STANDARD RESULT APPLIED / PROVED HERE for the accumulation.** Work at a fixed known warp. Every calibration and transport query lies in \([-D,D]\), \(D=R+3\), and all central queries lie in \([-2,2]\). On a fixed compact probability interval the empirical uniform CDF deviation \(\Delta_{a,n}\) obeys the interval Bernstein bound
\[
|\Delta(p)-\Delta(q)|
\le C\{\sqrt{|p-q|\log n/n}+\log n/n\}
\tag{9}
\]
with probability tending to one, uniformly in \(p,q\). Also
\(|\Delta(p)|\le C\{\sqrt{p(1-p)\log n/n}+\log n/n\}\);
maximum uniform order-statistic spacing is \(O_P(\log n/n)\).
These follow from binomial Bernstein bounds on rational grids and monotonicity, with a union bound and refinement. DKW supplies the compact \(O_P(n^{-1/2})\) supremum bound. These are CDF statements, not Hellinger bounds.

For the inverse CDF, the empirical overshoot is at most \(1/n\).
First localize its argument by DKW; apply (9) over the resulting
\(O_P(n^{-1/2})\) probability displacement and Taylor expand the true quantile. This gives, uniformly on the fixed query interval,
\[
\widetilde F_a(z)-F_a(z)
=-\frac{F'_a(z)}{\phi(z)}\Delta_{a,n}(\Phi z)
+O_{P,\delta,\theta}(N^{-3/4}\log N).
\tag{10}
\]
Polygon replacement uses the maximum-spacing bound and adds only
\(O_{P,\delta,\theta}(\log N/N)\).
The constants in this argument are finite public-envelope quantities, not unknown favorable radii: on \([-D-1,D+1]\) use
\[
q_a'=F_a'/\phi,\qquad
q_a''=(F_a''+zF_a')/\phi^2,
\tag{11}
\]
and bound \(F_a',F_a''\) from the original \(2e^x,10e^x\) envelopes and the explicit \(h',h''\). Thus the constants involve inverse powers of \(\phi(D+1)\). They are not bounded as \(\delta\downarrow0\).

For an explicit choice, write \(d=D+1\), \(H_d=1+\delta+2\omega d\), \(E_d=d+\omega d^2+a+1.1\), with \(\omega,a\) as in conditioning (1), and take
\[
F_{1,d}=2(1+H_d)e^{E_d},\quad
F_{2,d}=(10+10H_d^2+4\omega)e^{E_d},\quad
C_D=1+\frac{F_{1,d}}{\phi(d)}
+\frac{F_{2,d}+dF_{1,d}}{\phi(d)^2}.
\]
The compact inverse-Taylor/spacing remainder in (10) is bounded in probability by a universal multiple of \(C_D N^{-3/4}\log N\), once the compact localization event holds. This exposes the otherwise hidden calibration-tail constants.


A central difference uses at most \(12K\) signed response terms. Each transported evaluation has additionally \(2|m|+2\) normalizing/transport terms. Consequently its total nonlinear remainder is
\[
O_{P,\delta,\theta}((K+L+1)N^{-3/4}\log N)
=O_{P,\delta,\theta}(N^{-3/4}\log^2N),\qquad L=\lceil R+4\rceil.
\tag{12}
\]
Its population terminal bias is \(O_{\delta,\theta}(\lambda^K)=O_{\delta,\theta}(N^{-2})\). The kernel increment proof in conditioning (5) bounds its omitted infinite-series norm by \(D_c\lambda^{K/2}=O_{\delta,\theta}(N^{-1})\). Its omitted empirical average is smaller by \(N^{-1/2}\). These two separate controls justify growing \(K\); convergence of population telescoping alone would not.

For two coefficient arguments in one fixed chart at distance at most \(d\), apply (9) to paired queries. The six query maps are Lipschitz, and their \(k\)-th arguments move by at most \(2\lambda^k d\). Sum the square-root increments and smooth coefficients geometrically. Including finite transports, a uniform noise-increment bound is
\[
C_{\delta,\theta}\left\{
\frac{\sqrt{d\log N}}{\sqrt N(1-\sqrt\lambda)}
+\frac{(K+L+1)\log N}{N}
+\frac{d}{\sqrt N(1-\lambda)}
+(K+L+1)N^{-3/4}\log N\right\}.
\tag{13}
\]
Here \(C_{\delta,\theta}\) includes the finite compact constants and \(e^R\) row scaling. The finite number of transport terms is bounded by the same display after increasing that public constant. For \(d=B\log N/\sqrt N\) and fixed \(B,\delta,\theta\), (13) is \(o_P(N^{-1/2})\). Fixed chart membership is essential; global smoothness of an empirically reconstructed profile across arbitrary integer joins is not asserted.

## 4. All four coefficients and data-defined selection

**PROVED HERE.** The companion calibration proof gives global separation by \(\gamma\) and population contraction by \(q<1\), throughout the enlarged search box and full original profile class. The crude empirical calibration noise is
\(\varepsilon_N=O_{P,\delta,\theta}((K+L+1)/\sqrt N)\).
Projection is nonexpansive in infinity norm, so
\[
\|c_i^t-c_i\|_\infty
\le q^t\operatorname{diam}_\infty([.4,1.1]^2)
+\tau\varepsilon_N/(1-q).
\tag{14}
\]
After \(k_N\) steps the iterates lie in an \(O_{P,\delta,\theta}(\log N/\sqrt N)\) neighborhood. All original truths have distance at least .1 from the search boundary.

Only for analysis put \(\widehat R_i=\widehat y_i-\widehat B_i\),
\(J_i=\nabla B_i(c_i)\), and
\[
c_i^*=c_i+J_i^{-1}\widehat R_i(c_i).
\]
The fixed-query infinite-series expansion from Section 3 gives
\(c_i^*-c_i=O_P(N^{-1/2})\). The original \(g''\) envelope bounds the population second derivative on the fixed calibration domain. Taylor expansion and (13) therefore give
\(\widehat R_i(c_i^*)=o_P(N^{-1/2})\).
On a localized radius \(B\log N/\sqrt N\) event, let \(\xi_N\) denote the bound in (13). Nonexpansive projection, population contraction, and \(c_i^*\) lying in the box give
\[
\|c_i^{t+1}-c_i^*\|_\infty
\le q\|c_i^t-c_i^*\|_\infty+
\tau\{\xi_N+\|\widehat R_i(c_i^*)\|_\infty\}.
\tag{15}
\]
After the second \(k_N\) steps, and then letting the localization constant \(B\) increase,
\[
\widehat c_i-c_i=J_i^{-1}\widehat R_i(c_i)+o_P(N^{-1/2}).
\tag{16}
\]
This holds jointly for both pairs, with shared anchor observations. No empirical contraction, exact empirical root, independence among queries, cross-fitting or unobserved Jacobian in the algorithm is used.

## 5. Response tails, lattice and coefficient substitution

**PROVED HERE.** Write \(\omega=\delta\sqrt\theta<.01\). Conditioning (1)–(3) derive the moments, including the necessary condition \(2p\omega<1\).
For each actual source,
\[
|F_a|+|F_a'|+|F_a''|
\le C_{\delta,\theta}(1+|z|)^m e^{\omega z^2+|z|}.
\]
For \(U=A_c-W_c\) or \(W_c\), its first two coefficient derivatives and the needed first \(z\) derivative, away from \(-1\), obey
\[
|U(z)|+|U'(z)|
\le C_{\delta,\theta}(1+|z|)^m e^{-z^2/2+C|z|},
\tag{17}
\]
uniformly on the search box. To see this, \(h'\ge1\) gives
\(|h^{-1}(h(z)-b)-z|\le |b|\le2.2\); differentiating these known inverse-density expressions introduces only bounded inverse derivatives and polynomial factors. The outward sums preserve Gaussian tails by the elementary sum bound in conditioning Section 4. The same proof applies to the coefficient derivatives.

For completeness, the required integrated-quantile bridge is explicit. For a step empirical quantile let
\(w(p)=U(\Phi^{-1}p)/\phi(\Phi^{-1}p)\) on \(|\Phi^{-1}p|\le B_N\), zero elsewhere, and \(V(p)=\int_0^p w\).
Layer-cake integration for positive responses, separately for positive and negative weights, gives
\[
\int_{-B_N}^{B_N}U(z)[Q_n(\Phi z)-F(z)]\,dz
=-\int_\mathbb R F'(z)\{V(H_n(\Phi z))-V(\Phi z)\}\,dz.
\tag{18}
\]
The common total weight cancels the endpoint term. This identity is applied to step quantiles before polygon replacement.

Since \(n\Phi(-B_N)=N^{1/4+o(1)}\), the Bernstein event in (9) makes empirical and population tail probabilities relatively close. Probability segments crossing the retained support then correspond to \(|z|\le B_N+1\) with comparable local densities. On smooth pieces,
\[
|w'(p)|\le C(1+|z|)^m e^{z^2/2+C|z|}.
\]
The integrated Taylor remainder in (18) is at most
\[
\frac{C\log N}{N}\int_{|z|\le B_N+1}
(1+|z|)^m e^{\omega z^2+C|z|}\,dz
=N^{-1+1.5\omega+o(1)}.
\tag{19}
\]
At each jump \(z_0\in\{-B_N,-1,B_N\}\), the crossing strip has probability width \(O(d_0)\), where
\(d_0^2=C\Phi(z_0)(1-\Phi(z_0))\log N/N\).
The extra error is bounded by
\(CF'(z_0)|[w]_{z_0}|d_0^2/\phi(z_0)\).
Thus cutoff jumps have order (19), and the fixed join gives \(O_P(\log N/N)\). Adjacent uniform spacings bound polygon replacement by
\(C(\log N/N)F'(z)/\phi(z)\); weighted integration gives (19) again.

The population tail uses the product envelope
\(e^{-(1/2-\omega)z^2+C|z|}\).
The omitted linear influence uses the individual-source Bochner envelope
\(e^{-(1/4-\omega)z^2+C|z|}\), and its empirical average has the extra factor \(N^{-1/2}\). The simultaneous orders are therefore

| Error | Order at each fixed warp |
|---|---|
| Central accumulated nonlinear quantile error | \(O_P(N^{-3/4}\log^2N)\) |
| Central population truncation | \(O(N^{-2})\) |
| Central omitted influence norm | \(O(N^{-1})\) |
| Population response tails | \(N^{-3/4+1.5\omega+o(1)}\) |
| Tail area, cutoff and polygon nonlinear error | \(N^{-1+1.5\omega+o_P(1)}\) |
| Omitted tail linear empirical average | \(N^{-7/8+1.5\omega+o_P(1)}\) |
| Finite lattice error on capped data | \(N^{-10+o(1)}\) |

Here \(N^{o(1)}\) denotes a public-envelope expression \(C(\log N)^m e^{C\sqrt{\log N}}\) at the fixed warp, not an assumed unknown remainder. For the last row, outward omitted arguments have magnitude at least \(M_N-O(1)\); their weights are \(N^{-12+o(1)}\). Capped responses and the \(O(K)\) finite reconstruction supply at most \(N^{2+o(1)}\). All rows needed for the statistic are \(o_P(N^{-1/2})\). In particular the old .1-warp tail exponents have been rederived, not reused.

Coefficient substitution also has a controlled remainder. The first two coefficient derivatives satisfy (17); the integrated absolute linear empirical noise is \(O_P(N^{-1/2})\), by
\(E|\Delta(p)|\le\sqrt{p(1-p)/n}\) and the integrable Bochner envelope. The central part has the same bound by geometric summability. Multiplication by (16) gives \(O_P(N^{-1})\). The population target's second derivative is bounded by the original profile envelopes and moments, giving an \(O_P(N^{-1})\) Taylor remainder. The finite truncation bounds apply to the derivative weights as well. Thus every first-order coefficient correction is retained.

## 6. Actual influence and full tangent closure

**PROVED HERE.** Use the source-vector kernels \(\kappa_{a,z}\) in conditioning (4), including their allocation factor \(1/\pi_a=5\). Apply the observable six-term functional to this vector:
\[
e_x=E(\kappa;x),\qquad
d_{x,r}=-\sum_{k\ge0}(e_{C^kx}-e_{C^kr}).
\]
For the fixed calibration charts set
\[
U_x=\kappa_{0,r}-\kappa_{3,r}+d_{x+m,r}+H_m(\kappa;x),
\qquad V_y=\kappa_{0,h^{-1}y}-U_{h^{-1}y}.
\]
Let \(z_-=-R,z_+=R,w_-=e^R,w_+=e^{-h(R)}\), and define
\[
\zeta_{i,\sigma}
=w_\sigma\{\kappa_{i,z_\sigma}
-U_{z_\sigma+\alpha_i}-V_{h(z_\sigma)+\beta_i}\},\qquad
\chi_i=J_i^{-1}(\zeta_{i,-},\zeta_{i,+})^\top .
\tag{20}
\]
The complete influence is
\[
\psi_A=\int_\mathbb R\{(A_c-W_c)\kappa_{0,z}
+W_c\kappa_{3,z}\}\,dz+\int_{-1}^0P_c(t)d_{t,r}\,dt,
\]
\[
\boxed{\psi_\eta=\psi_A+\sum_{i=1}^2 b_i^\top\chi_i,\qquad
b_i=(\partial_{\alpha_i}I,\partial_{\beta_i}I).}
\tag{21}
\]
The true quantities in (20)–(21) are analysis objects only. Sections 3–5 derive exactly this influence for (3), (7), rather than merely finding some solution of an adjoint equation.

Conditioning proves convergence of each infinite series and integral in the direct sum of individual source \(L^2\) spaces, before combining sources. In particular
\[
V_\eta=\sum_{a=0}^4\pi_a E_a\psi_{\eta,a}^2
\le\Psi(\delta,\theta)^2<\infty.
\]
All terms from source zero, including both calibration pairs and the target, are added before squaring. There are no independent copies of the reference quantile. The \(1/N\) normalization in (1) follows from \(5/N=1/n\).

For an actual admissible fixed profile/coefficient path write its curve velocity as \(q_a(z)\). Its density score in latent coordinates is
\[
S_a(z)=z\,u_a(z)-u_a'(z),\qquad u_a=q_a/F_a'.
\tag{22}
\]
Here \(q_a\) includes both analytic profile directions and every coefficient derivative for sources 1 and 2. Differentiating the monotone pushforward density proves (22). Boundary terms vanish by the Gaussian factor and the envelopes. Pairing a quantile kernel with this score gives
\(\langle\kappa_{a,z},S\rangle_\pi=q_a(z)\): integrate
\((-\phi u_a)'=\phi S_a\) from \(-\infty\) to \(z\).
The commutator identity then recovers the profile-direction differences; its normalizer recovers \(r_1(r)\) because \(r_1(0)=0\).
For a pure profile direction, (20) is zero. For coefficient directions it equals \(J_i\,d_i\) before multiplication by \(J_i^{-1}\). Equation (6) therefore gives
\[
\langle\psi_\eta,S\rangle_\pi=\dot I_\eta
\]
for all such paths. A finite dictionary has not been substituted for this class.

At profile-envelope boundary truths, not every formal affine perturbation remains in the model. The path domain consists of actual two-sided admissible weighted-\(C^3\) paths, with derivatives \(r_j(0)=0\), \(|r_j^{(k)}(x)|\le C_v e^x\), \(k\le3\), and finite coefficient directions; their closure in score norm is the tangent closure asserted here. At interior truths these include the usual allowed analytic affine directions. The continuous pairing with the finite-norm \(\psi_\eta\) extends the derivative to that full closure. No claim of additional directions at a constrained boundary is needed.

**STANDARD RESULT APPLIED, hypotheses checked.** These fixed smooth admissible paths are DQM: the lower derivative envelope bounds \(q_a/F_a'\) and its derivatives by polynomials in \(|z|\); the upper profile derivative and warp envelopes control the differentiated pushforward square-root density in Gaussian \(L^2\). Taylor differentiation on compact \(z\)-sets followed by Gaussian dominated tail control gives the DQM remainder. For a nonlinear path this uses differentiability in the stated weighted \(C^3\) norm; an arbitrary pointwise differentiable path would be insufficient. The fixed five-group product experiment is LAN. The joint CLT of the five-source score and (21), followed by change of measure, transfers the \(o_P(N^{-1/2})\) remainder along each fixed contiguous \(t/\sqrt N\) path and cancels the target's local mean shift. This proves fixed-path regularity. It implies neither uniformity over changing directions nor a neighborhood-uniform CLT.

More explicitly, \(q_a/F_a'\) is uniformly bounded and its first two derivatives have polynomial envelopes. Nearby admissible paths shift the inverse latent coordinate by \(O(|t|)\) uniformly, by integrating the bounded velocity-to-derivative ratio. The differentiated square-root densities are consequently dominated by a polynomial in \(|z|\) times \(e^{C|z|}\sqrt{\phi(z)}\), an \(L^2(dz)\) envelope. This justifies the stated compact-to-tail passage.


## 7. A finite certified arithmetic specification

**PROVED HERE as a mathematical realization, not an implemented backend.** A simple midpoint rule suffices; the old fast knot-integration or bit-complexity theorem is not imported. Inputs are observed responses and known warp constants supplied through refinable certified enclosures. More working precision cannot recover digits absent from a fixed-resolution observation scheme.

The following intentionally loose public bounds give explicit tolerances. Set \(B=B_N,M=M_N,Z_*=B+M+3\),
\[
H=1+\delta+2\omega Z_*,\quad
\ell_0=4(1+H),\quad \ell_1=8(1+H^2),\quad
a_0=4H,\quad a_1=4+8H^2,
\]
\[
U_0=a_0+2M\ell_0,\quad U_1=a_1+2M\ell_1,\quad
P_0=(2M+1)\ell_0,\quad P_1=(2M+1)\ell_1.
\tag{23}
\]
Indeed \(|\phi|,|\phi'|\le1\), \(\rho_2\le1\), \(|\rho_2'|\le2\), \(h'\ge1\), and \(|h''|<1\); differentiating \(A,\mathcal L\) gives these bounds on all finite queries. Thus \(U_0,U_1\) bound the sum of absolute outer weights and their derivatives, and \(P_0,P_1\) bound the finite central weights and their derivatives.

On positive capped data the polygon is at most \(N^2\) and \(N^3\)-Lipschitz in probability. The finite central difference is at most \(12KN^2\) and \(12KN^3\phi(0)\)-Lipschitz in \(t\), since every varying central query is 2-Lipschitz. Hence the outer and central integrands are Lipschitz, separately across the known join, with bounds
\[
L_o=N^2U_1+N^3\phi(0)U_0,\qquad
L_c=12KN^2P_1+12KN^3\phi(0)P_0.
\]
Split the outer integral at \(-1\) and use midpoint panels of length at most
\[
h_N^{\rm mesh}=[16N(1+2B L_o+L_c)]^{-1}.
\tag{24}
\]
The absolute midpoint error is at most one half the mesh times the sum of interval lengths times their Lipschitz bounds, hence at most \(1/(32N)\). Polygon knots create no jump and require no derivative of an unknown profile.

For height sensitivity put
\[
A_{\rm cal}=6+24K+4L,\quad
C_{\rm cal}=\tau e^R A_{\rm cal},\quad
J_Y=2B U_0+12KP_0,
\]
and read each height to absolute error
\[
\epsilon_{Y,N}\le[64N(1+C_{\rm cal}+J_Y)]^{-1}.
\tag{25}
\]
There are at most \(A_{\rm cal}\) signed quantile evaluations in a calibration residual row, including its two normalizations and transports. Sorting is 1-Lipschitz in maximum-coordinate distance; thus the polygon changes by at most \(\epsilon_{Y,N}\), even if order changes. Equation (25) contributes at most \(1/(64N)\) to each update and to the fixed-coefficient target.

Require each complete numerical update, at its current representable coefficient argument, to be enclosed to error at most \(1/(4N)\) in infinity norm relative to the exact-data update at that same argument. Include input error, known row scalings, \(\tau\), finite operators, projection, starting-state and summation rounding. For example four shares of \(1/(16N)\) suffice. For the target require the sum of integration, input, primitive, node/endpoint and final summation errors to be at most \(1/(4N)\), with (24) and (25) using less than their shares.

After rounding, clamp each stored coefficient to the exact rational box \([.4,1.1]^2\); the executed states therefore remain in the box. Endpoint/node enclosures and final clipping roundoff belong to the displayed target share, rather than being uncharged extra operations.


These are attainable finite operations: inverse-warp bisection uses \(h'\ge1\), square root/exponential/Gaussian probabilities admit certified rational approximations, and every finite expression is continuous on its compact coefficient domain. Query errors are converted to response errors using \(N^3\); finite sums admit absolute-error allocation term by term. Projection and polygon evaluation are continuous Lipschitz maps, including at knots. An interval crossing a knot is evaluated over the whole interval rather than requiring an undecidable exact rank comparison. Certified refinement terminates for any positive absolute tolerance. Public chart preprocessing was made finite above. For schedule ceilings use outward enclosures and choose a modestly larger integer with certified \(\lambda^K<N^{-2}\), \(q^{k_N}<N^{-1}\); the added fixed slack changes none of the proof. This defines a finite statistic without an oracle accuracy test on an unknown truth.

This proof does not retain a claimed \(N\,\mathrm{polylog}N\) backend bound or a universal small working-bit constant. Constants include the actual extreme calibration queries. No midpoint panels or certified primitives were executed.

## 8. Statistical arithmetic comparison and conclusion

The numerical errors in Section 7 do not require empirical contraction. Add \(1/(4N)\) to (14) and (15). They preserve localization and the same expansion (16), so the computed and ideal four-vectors differ by \(o_P(N^{-1/2})\). This is an in-probability conclusion, not deterministic \(O(N^{-1})\) agreement of their whole trajectories.

The finite target's coefficient sensitivity is \(O_P(1)\), as follows directly. The central empirical reconstruction is its bounded population value plus \(O_P(K/\sqrt N)\), and the coefficient-derivative periodized weights have uniformly bounded integrals. For the outer part positivity and an order-statistic count give
\[
E\widetilde Q_a(p)\le M_4^{1/4}(1-p)^{-1/4},\qquad 0\le p<1.
\tag{26}
\]
To check this for the polygon, bound it by the order statistic at rank \(\max(1,\lceil np\rceil)\); at least \(n(1-p)\) observations are no smaller than that order statistic. Apply the fourth-moment bound and Jensen. The coefficient-derivative version of (17), multiplied by \((1-\Phi z)^{-1/4}\), is integrable: its positive-tail quadratic exponent is at most \(-3z^2/8+C|z|\). Tonelli and Markov prove the claimed tight gradient bound, uniformly on the public coefficient box.

Decompose the computed target minus \(H_N(\widetilde F,\widehat c)\) into numerical integration error, direct height perturbation at the computed coefficients, and coefficient perturbation with exact heights. Sections 7 and (26) bound these by \(O(N^{-1})\), \(J_Y\epsilon_{Y,N}=O(N^{-1})\), and \(O_P(1)o_P(N^{-1/2})\), respectively. Thus the complete arithmetic statistic is asymptotically equivalent to the ideal one.

For a finite-bit gate check encoded heights for positivity and the rational cap \(N^2\). Both exact and encoded gates accept whenever
\(\epsilon_{Y,N}<\min Y\) and \(\max Y<N^2-\epsilon_{Y,N}\).
Each source has \(Y\ge\tfrac12e^Z\), so
\[
P(\min Y\le\epsilon_{Y,N})\le N\Phi(\log(2\epsilon_{Y,N})),\qquad
P(\max Y\ge N^2-\epsilon_{Y,N})
\le NM_4/(N^2-\epsilon_{Y,N})^4.
\]
Both vanish; the first does so faster than any fixed inverse power at the fixed warp. Hence no exact finite-bit classification of an arbitrary real at the cap is assumed. Clipping at (4) is 1-Lipschitz and asymptotically inactive. These negligible events also transfer along the already specified fixed contiguous paths.

Combining Sections 3–8 proves (1) for a feasible source-data statistic in the explicit arithmetic input model. The coefficient roots, profile derivatives, scores and influence are used only to analyze it.

**OPEN:** sharp separation rates, shrinking-\(\delta\) uniformity, finite-sample MSE or accuracy, efficient variance, feasible variance estimation and studentization, practical onset, fixed-resolution-data correction and a certified executable backend. In particular a bounded clipped estimator and a pointwise CLT alone do not prove \(O(1/N)\) MSE. The historical backend no-go and unresolved efficiency decisions are unchanged.
