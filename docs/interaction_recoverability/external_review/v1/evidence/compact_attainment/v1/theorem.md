# Pointwise attainment and computational comparison

**PROVED HERE, in the explicit arithmetic model of [quadrature_proof.md](quadrature_proof.md):** the compact statistic admits a pointwise regular root-\(N\) expansion. Public projected calibration and piecewise certified integration give a mathematical realization with \(N\,\mathrm{polylog}N\) bit operations and \(O((\log N)^2)\) precision. The supplied diagnostic code does not implement the complete interval backend. The result does not establish affordable constants, a variance improvement, efficiency or MSE.

## The ideal and accelerated statistics

Let \(\widetilde Q_a\) be the continuous polygon through the observed order statistics at \(i/n\), constant on \([0,1/n]\), and \(\widetilde F_a(z)=\widetilde Q_a(\Phi z)\). The construction uses these five groups only. Select
\[
K=\lceil 2\log N/|\log(271/290)|\rceil+2,\quad
R^2=\tfrac32\log N,\quad M=\lceil4R+30\rceil+2. \tag{1}
\]
At \(N=5n\), \(n\ge2\), define the **compact ideal** by exact integration of construction (6), with the lexicographically first minimum-residual coefficient pair on each public rectangle grid of mesh at most \(N^{-1}\). This proof comparator is finite and data-defined; it is not executed.

Define the **compact iteration statistic** by replacing those two selections with the prescribed public iterations below. Define its arithmetic realization by the explicit rounding, knot enclosure, gap exclusion, two-node Gauss and arithmetic allocations in the companion proof. All three use the same \(K,M,R\), observed response-cap fallback and clipping to \([-1000,1000]\). The cap event is checked on the original observable responses. Inside it, upward dyadic rounding clipped at \(N^2\) preserves positivity and the cap; this avoids assuming that an independently rounded gate makes identical decisions. If input precision is externally fixed, this theorem's rounding comparison does not automatically apply.

The conclusions at each fixed truth in the declared local class are
\[
\widehat I_{\mathrm{compact\ ideal}}-I_\eta
=\frac1N\sum_{a=0}^4\sum_{i=1}^n\psi^{\,\mathrm{new}}_{\eta,a}(Y_{ai})
+o_{P_\eta}(N^{-1/2}), \tag{2}
\]
\[
\sqrt N(\widehat I_{\mathrm{iteration}}-\widehat I_{\mathrm{compact\ ideal}})
\ \longrightarrow_P0,\qquad
|\widehat I_{\mathrm{fast}}-\widehat I_{\mathrm{iteration}}|
=O(N^{-1})\quad\text{on the cap event}. \tag{3}
\]
The first comparison is statistical; the second is a deterministic arithmetic comparison. Neither comparator is the old anchor statistic.

## Compact empirical-process control

**STANDARD RESULT APPLIED, with the growing-series argument checked here.** On the fixed compact response-query interval \(D=[-7/4,13/4]\), the original proof's uniform empirical interval bound and inverse-CDF expansion give
\[
\widetilde F_a(z)-F_a(z)
=-\frac{F'_a(z)}{\phi(z)}\Delta_{a,n}(\Phi z)
+O_P(N^{-3/4}\log N), \tag{4}
\]
uniformly in \(z\). The probability-integral transform \(\Delta_{a,n}\) is an analysis variable, never an estimator input. Polygon interpolation changes step quantiles by \(O_P(\log N/N)\) on this fixed compact set. DKW/interval concentration controls CDFs and quantile remainders here; it is not asserted to bound Hellinger distance.

There are at most \(36K\) responses in a central difference, so the total nonlinear remainder is \(O_P(N^{-3/4}\log^2N)\). Population truncation in construction (3) is \(O(N^{-2})\). For each signed response query, its centered quantile-kernel increment has \(L^2\) norm bounded by a public constant times the square root of the distance. Contraction makes the infinite central series converge geometrically in each individual source \(L^2\). The omitted series norm is \(O(\lambda^{K/2})=O(N^{-1})\). These two facts justify growing \(K\); a fixed-\(K\) CLT alone would not.

For arguments in any one verified calibration chart with distance at most \(d_N\), summing the empirical interval increments geometrically gives the noise increment
\[
O_P\{N^{-1/2}\sqrt{d_N\log N}+K\log N/N+d_N/\sqrt N\}
+O_P(N^{-3/4}\log^2N). \tag{5}
\]
This is \(o_P(N^{-1/2})\) for \(d_N=B\log N/\sqrt N\).
The finitely many normalization/transport terms have the same or better modulus. Chart membership was proved over the full public rectangles in construction.md. No global continuity across integer joins is required or asserted.

## Public projected calibration

For pair \(i\), let \(M_{i0}\) be the public matrix in the request and \(A_i=M_{i0}^{-1}\). The true population map has Jacobian
\[
J_i(c)=M_{i0}\operatorname{diag}(e^{c-c_{i0}})+E_i(c),\qquad
|E_{i,uv}|\le.01 e^{.02}(M_{i0})_{uv}.
\]
It follows, using operator norm bounded by Frobenius norm for the error, that
\[
\|I-A_iJ_i(c)\|_2\le q_i=
e^{.02}-1+.01e^{.02}\|A_i\|_F\|M_{i0}\|_F. \tag{6}
\]
The relevant Frobenius product equals
\[
\frac{e^{\alpha_{i0}-\beta_{i0}}\cosh 2+
e^{\beta_{i0}-\alpha_{i0}}\cosh(2h_1)}
{\sinh(h_1-1)}.
\]
Since \(|\alpha_{i0}-\beta_{i0}|\le1/15\), \(h_1<1.15\), and \(h_1-1>.14\),
this is less than
\((15/14)(89/10)/(7/50)=6675/98<69\).
Here \(e^{1/15}<15/14\), \(\cosh2<19/5\), and \(\cosh(23/10)<51/10\).
For a purely rational certificate of the latter two bounds, sum the cosine-hyperbolic series through degree 10 and bound the first omitted term and following terms by a geometric series with ratio at most \( (23/10)^2/(13\cdot14)<.03\).
Using \(e^{.02}<50/49\) yields
\[
q_i<1/49+69/98=71/98<\bar q=3/4. \tag{7}
\]
Convex-rectangle integration of the derivative in (6), then nonexpansiveness of projection, proves contraction of the **population** update.

Set \(k_N=\lceil\log N/|\log(3/4)|\rceil\). Start at the public center and execute \(2k_N\) updates
\[
c^{t+1}=\Pi_{\mathcal R_i}\{c^t+A_i[\widehat y_i-\widehat B_i(c^t)]\}.
\]
The crude uniform response/reconstruction error is
\(\varepsilon_N=O_P(\log N/\sqrt N)\). Thus
\[
\|c^t-c_i\|\le \bar q^t\operatorname{diam}(\mathcal R_i)
+\|A_i\|\varepsilon_N/(1-\bar q). \tag{8}
\]
After \(k_N\) steps the iterates are in the neighborhood required for (5).

Only in the proof, put
\(c_i^*=c_i+J_i(c_i)^{-1}\widehat R_i(c_i)\),
\(\widehat R_i=\widehat y_i-\widehat B_i\).
The fixed-query expansion shows \(c_i^*-c_i=O_P(N^{-1/2})\).
Taylor expansion of the population map and (5) give
\(\widehat R_i(c_i^*)=o_P(N^{-1/2})\).
On the localized neighborhood, the supremum of empirical-noise increments is \(\xi_N=o_P(N^{-1/2})\). Nonexpansive projection and population contraction then give
\[
\|c^{t+1}-c_i^*\|\le\bar q\|c^t-c_i^*\|
+\|A_i\|[\xi_N+\|\widehat R_i(c_i^*)\|].
\]
A second \(k_N\) steps prove
\[
\widehat c_i-c_i=J_i(c_i)^{-1}\widehat R_i(c_i)+o_P(N^{-1/2}). \tag{9}
\]
The compact ideal grid has the same expansion: the population separation argument from the saved theorem applies to the unchanged \(B_i\), the nearest grid point to \(c_i^*\) has residual \(o_P(N^{-1/2})\), and (5) controls the actual selected residual. This proves coefficient equivalence. Neither the true root nor \(J_i(c_i)\) is an algorithm input, and the empirical update was never assumed contractive. Uniformly bounded numerical update errors \(o(N^{-1/2})\) can also enter these recurrences; the companion proof instead supplies a stronger deterministic comparison of the complete finite iteration.

## Actual new influence, all sources and coefficients

For a source label \(a\), define the centered stratified vector whose only nonzero coordinate is
\[
\kappa_{a,z}(y)=5\frac{F_a'(z)}{\phi(z)}
  \{\Phi(z)-{\bf1}[y\le F_a(z)]\}. \tag{10}
\]
Let \(e_x=E_{45}(\kappa;x)\), taking the signed source-specific evaluation pattern in construction (1), and
\[
d_{x,r}=-\sum_{k\ge0}(e_{C^kx}-e_{C^kr}),\qquad
u_x=\kappa_{0,r}-\kappa_{3,r}+d_{x+m(x),r}+H_{m(x)}(\kappa;x),
\]
\[
v_y=\kappa_{0,h^{-1}y}-u_{h^{-1}y}.
\]
The infinite series converges in each source space by construction (7). Define the residual vector and the pair coefficient duals
\[
\zeta_i=
\begin{pmatrix}
\kappa_{i,-1}-u_{-1+\alpha_i}-v_{-h_1+\beta_i}\\
\kappa_{i,1}-u_{1+\alpha_i}-v_{h_1+\beta_i}
\end{pmatrix},\qquad
\chi_i=J_i(c_i)^{-1}\zeta_i. \tag{11}
\]
The anchor target influence is
\[
\psi_A^{\rm new}=
\int_{\mathbb R}\{(A-W_b)\kappa_{0,z}+W_b\kappa_{3,z}\}\,dz
+\int_J P(t)d_{t,r}\,dt. \tag{12}
\]
Writing \(b_i=(\partial_{\alpha_i}I,\partial_{\beta_i}I)\) in pair order, the complete influence is
\[
\boxed{\ \psi_\eta^{\rm new}=\psi_A^{\rm new}+\sum_{i=1}^2 b_i^\top\chi_i.\ } \tag{13}
\]
The derivatives in \(b_i\) are those of the actual signed target with the profiles fixed. For example \(\partial_{\alpha_1}K_1=-\phi'(\cdot-\alpha_1-\alpha_2)+\phi'(\cdot-\alpha_1)\). Profiles, derivatives, Jacobians and influence values in (10)–(13) are proof objects; computing the estimator does not require them.

Each of sources \(0,1,2,3,4\) appears in (13). The asymptotic variance is
\[
V_\eta^{\rm new}=\sum_{a=0}^4\tfrac15 E_\eta[(\psi_{\eta,a}^{\rm new})^2].
\]
All repeated contributions from a given group, especially source 0, are added before squaring. There is no independence assumption among different quantile queries from a group and no hidden cross-fitting.

The profile and coefficient score equations follow by pairing (10) with the relevant five-law tangent scores and applying the population/tangent identities. Thus (13) represents the full target derivative on the closure of the admissible tangent space, not merely on a finite dictionary. Its off-model extension depends on the new anchor and unit cell; equality with the old influence or variance is unproved and is not used.

## Tails, lattice and nonlinear terms

The new weights and their first two coefficient derivatives obey the same type of integrable envelopes as the saved proof:
\[
|U(z)|+|U'(z)|\le C(1+|z|)^m e^{-z^2/2+C|z|},
\quad F_a'(z)\le C(1+|z|)^m e^{.1z^2+C|z|}.
\]
The new fixed jump is at \(b=-7/4\). Moving the boundary alters finite constants, not Gaussian tail exponents. In particular the sourcewise Bochner envelope for \(U(z)\kappa_{a,z}\) is \(C(1+|z|)^m e^{-.15z^2+C|z|}\), so every individual source \(L^2\) norm is finite before sources are combined.

Apply the empirical-CDF area identity of the saved theorem separately to the positive and negative signed weights for step quantiles, then bound their polygon replacement by the adjacent-order-statistic spacing. There are three jumps: \(-R,b,R\). The fixed jump contributes \(O_P(\log N/N)\); the cutoff jumps satisfy the same tail bound as the smooth remainder. The resulting simultaneous orders are:

| Term | Bound |
|---|---|
| Compact nonlinear quantile remainder | \(O_P(N^{-3/4}\log^2N)\) |
| Compact population series truncation | \(O(N^{-2})\) |
| Compact omitted influence norm | \(O(N^{-1})\) |
| Population response tails | \(N^{-.6+o(1)}\) |
| Tail area/interpolation nonlinear error | \(N^{-.85+o_P(1)}\) |
| Omitted tail linear empirical term | \(N^{-.725+o_P(1)}\) |
| Finite lattice error under response cap | \(N^{-10+o(1)}\) |

For clarity, the tail calculation is \((\log N/N)\int_{|z|\le R+1}e^{.1z^2+C|z|}(1+|z|)^m dz=N^{-.85+o(1)}\); bias uses \(\int_{|z|>R}e^{-.4z^2+C|z|}(1+|z|)^m dz=N^{-.6+o(1)}\). With \(R^2=c\log N\), the required inequalities \(c<2\), \(.4c>.5\), and \(.1c<.5\) all hold at \(c=1.5\).
Here \(N^{o(1)}\) means a fixed public-envelope expression \(C(\log N)^m e^{C\sqrt{\log N}}\), not an unknown favorable radius. The outward direction makes each omitted lattice argument have magnitude at least \(M-5\), even when \(R\) grows. Shifted Gaussian densities then give \(N^{-12+o(1)}\) omitted weights; multiplication by the response cap and \(K=O(\log N)\) gives the stated \(N^{-10+o(1)}\) bound; the additive 30 affects constants, not a claimed finite-sample threshold.

Coefficient substitution is also controlled. Uniform first and second coefficient-derivative weight envelopes are integrable. The absolute weighted empirical linear noise times \(\|\widehat c-c\|\) is \(O_P(N^{-1})\), and the population Taylor remainder is \(O_P(N^{-1})\). Therefore (9) adds exactly the four terms in (13), without unbounded nonlinear remainders.

The unchanged response fourth-moment bound gives
\(P(\max Y>N^2)\le M_4N^{-7}\), where
\(M_4=128e^4(e^8+\sqrt5e^{40.2})\).
The target magnitude is below 160 in the original envelope, so clipping at 1000 is asymptotically inactive. These are asymptotic fallback justifications; their large constants are not practical risk bounds.

## Arithmetic comparison, local alternatives and limits

The [arithmetic proof](quadrature_proof.md) bounds fourth derivatives of smooth polygon pieces using only known warp/density derivatives, public coefficient rectangles and the response cap. It encloses all knot preimages, excludes tiny gaps with an absolute integral bound, splits at \(b\), and allocates Gauss, primitive, summation, input and coefficient errors explicitly. It proves \(O(N^{-1})\) approximation to the finite-iteration statistic on capped data with input and working precision \(O((\log N)^2)\) and total bit work \(N\,\mathrm{polylog}N\). This precision is a sufficient worst-case bound, not a necessary lower bound. Real observations must be read/rounded to the specified increasing accuracy; a fixed machine-precision recording scheme is not certified by this argument.

Under each fixed admissible DQM path \(r_j(0)=0\), \(|r_j^{(k)}|\le C_v e^x\), \(k\le3\), and four finite coefficient directions, the saved five-law LAN argument applies. The joint CLT for (13) and the full score, followed by change of measure, transfers the remainder by contiguity and cancels the target's local mean shift. This gives fixed-path \(N^{-1/2}\) regularity at truths for which the path remains admissible. It does not imply uniformity over changing directions or all nearby truths.

**OPEN:** efficient variance attainment, feasible variance estimation/studentization, \(O(1/N)\) MSE, neighborhood-uniform coverage, global minimax rates, optimal computational complexity, and affordable validated performance. No CLT is promoted to any of these. The mathematical arithmetic realization has been specified; certification of the full delivered backend remains implementation work.

