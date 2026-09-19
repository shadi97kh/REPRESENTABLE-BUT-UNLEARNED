# Local attainment with five observed laws — v1

**Decision: local attainment established for the mathematical estimator specified below.** The finite numerical realization is analyzed under an explicit certified-arithmetic model. The supplied, sample-disabled arbitrary-precision implementation is a reference implementation, not a certified interval backend. No sampled data were generated or processed in this continuation.

## Claim and experiment

**PROVED HERE.** Retain the original [model](../../unknown_profiles/model.md), known Gaussian latent law, known \(h(z)=z+.1z\sqrt{1+z^2}\), and the five independent source groups \(0,e_1,e_2,e_3,e_4\), with \(n\) observations each and \(N=5n\). The latent variables are unobserved. Let the public reference centers be \(\alpha_0=(2/3,3/4)\), \(\beta_0=(3/5,4/5)\).

The theorem's truth set \(\mathcal U\) consists of original-class profiles satisfying
\[
\max_{j=1,2;\,0\le k\le3}\sup_x e^{-x}|g_j^{(k)}(x)-e^x|\le .01,
\qquad \|c-c_0\|_\infty\le .01 .
\tag{1}
\]
This is a declared parameter neighborhood. The statistic receives only the five response lists and their sizes. It receives no true profile, coefficient, response curve, derivative, latent value, target, or influence function. Public coefficient search rectangles have radius .02 around the reference centers.

For every fixed truth \(\eta\in\mathcal U\), the one statistic in [estimator.md](estimator.md) satisfies
\[
\widehat I_N-I(\eta)=\frac1N\sum_{a=0}^4\sum_{i=1}^n
 \psi_{\eta,a}(Y_{ai})+o_{P_\eta}(N^{-1/2}).
\tag{2}
\]
The centered source influences are individually square integrable and are precisely the explicit, generally nonminimum-norm solution (25) in the [local proof](../../local_five_law/v1/proof.md), evaluated at \(\eta\). Thus
\[
\sqrt N(\widehat I_N-I(\eta))\ \Rightarrow\ N(0,V_\eta),\qquad
V_\eta=\tfrac15\sum_{a=0}^4E_\eta\psi_{\eta,a}^2.
\tag{3}
\]
At the reference this estimator is regular along every stipulated fixed admissible analytic \(N^{-1/2}\) path, with both profiles and all four coefficients perturbed. The same holds at truths with a margin within the declared neighborhood. These are pointwise CLTs and fixed-path local regularity, not a asserted uniform CLT over all sequences in \(\mathcal U\), MSE theorem, efficient bound, or global minimax theorem.

## Population identity and the finite statistic

Use the notation and checked kernels of the local proof, now with **population response curves** \(F_a\), not its tangent \(q_a=\dot F_a\). Define \(T_jx=h^{-1}(h(x)+j)\), \(Cx=T_2x-1\), and
\[
\mathcal E_F(x)=F_0(T_2x)+F_0(T_1x)+F_0(Cx)
 -F_4(x)-F_4(T_1x)-F_3(Cx)=g_1(Cx)-g_1(x).
\]
On \(J=[3,6]\), \(C(J)\subset J\), \(0<C'\le\lambda=575/613<1\).
Put
\[
D_K(F;x,y)=-\sum_{k=0}^{K-1}
 [\mathcal E_F(C^kx)-\mathcal E_F(C^ky)].
\tag{4}
\]
Then \(D_K=g_1(x)-g_1(y)-[g_1(C^Kx)-g_1(C^Ky)]\). The last term is bounded by \(2e^6\lambda^K|x-y|\). Normalization gives
\[
g_1(4)=1+\sum_{j=0}^3[F_3(j)-F_0(j)].
\tag{5}
\]
Equation (5), (4), and finite unit translations give a finite \(\widehat g_{1,K}(x)\); set \(\widehat g_{2,K}(y)=\widetilde F_0(h^{-1}y)-\widehat g_{1,K}(h^{-1}y)\). The hat denotes the same algebra applied to the data curves \(\widetilde F\). No analyticity of the empirical reconstruction is required.

For any candidate coefficients \(c\), write \(A_c=h'K_{2,c}\), \(L_c=K_{1,c}-A_c\),
\[
P_c(t)=\sum_{j\in\mathbb Z}L_c(t+j),\quad
W_c(z)=
\begin{cases}\sum_{j\ge1}L_c(z+j),&z\ge4,\\
-\sum_{j\ge0}L_c(z-j),&z<4.
\end{cases}
\]
The original target at \((g,c)\) equals
\[
\int[(A_c-W_c)F_0+W_cF_3]\,dz+
\int_4^5P_c(t)[g_1(t)-g_1(4)]\,dt .
\tag{6}
\]
It is important that (6) holds for every candidate \(c\), while its anchor response curves do not depend on \(c\).

The estimator uses polygonal empirical quantiles \(\widetilde Q_a\), specified below, \(\widetilde F_a(z)=\widetilde Q_a(\Phi z)\), two minimum-residual coefficient selections, and (6) with cutoff \([-R_N,R_N]\), finite lattices \(M_N\), and \(D_{K_N}\) in the central integral. It clips the result to \([-1000,1000]\).

## A. Empirical quantiles and a growing compact reconstruction

**STANDARD RESULT APPLIED / elementary proof supplied.** In analysis only, transform observations to independent uniforms \(U_{ai}=P_{\eta,a}(Y_{ai})\), with empirical CDF \(H_{a,n}\). They are not estimator inputs. Put \(\Delta_{a,n}(p)=H_{a,n}(p)-p\).

For every fixed desired polynomial failure exponent, binomial Bernstein and a union bound give
\[
|\Delta_{a,n}(p)|\le C\{\sqrt{p(1-p)\log n/n}+\log n/n\},
\tag{7}
\]
and, simultaneously for every probability interval \(I\),
\[
|(P_n-P)I|\le C\{\sqrt{|I|\log n/n}+\log n/n\}.
\tag{8}
\]
For a direct justification, take a grid of mesh \(n^{-3}\), use Bernstein for its at most \(O(n^6)\) endpoint pairs, and enclose arbitrary intervals between neighboring grid intervals. Endpoint strips satisfy the same bound. Increase the fixed logarithmic constant to obtain any prescribed \(O(n^{-A})\) failure probability. The same argument with empty intervals gives maximum uniform spacing \(O(\log n/n)\). This grid is only a proof device.

Let \(Q_n(p)=Y_{(\lceil np\rceil)}\) for \(0<p\le1\), with \(Q_n(0)=Y_{(1)}\). If \(v=H_n^{-1}(p)\), the rank overshoot is at most \(1/n\), and
\[
v-p=-\Delta_n(p)-[\Delta_n(v)-\Delta_n(p)]+O(1/n).
\]
On a fixed interior probability interval, DKW bounds \(|v-p|=O_P(\sqrt{\log n/n})\), and (8) bounds the bracket by \(O_P(n^{-3/4}\log n)\). Ordinary \(C^2\) Taylor expansion of \(F(\Phi^{-1}p)\) gives, uniformly on every fixed compact latent interval needed by the construction,
\[
Q_n(\Phi z)-F(z)=
-\frac{F'(z)}{\phi(z)}\Delta_n(\Phi z)
+O_P(n^{-3/4}\log n).
\tag{9}
\]
All response derivatives and inverse-density bounds on such compact intervals are bounded by known class/warp envelopes. No empirical profile derivative is used.

The implemented quantile convention is the continuous polygon through \((i/n,Y_{(i)})\), \(i=1,\ldots,n\), constant below \(1/n\) and above 1. It differs from \(Q_n(p)\) by at most one adjacent order-statistic spacing. The uniform-spacing event and the mean value theorem give \(O_P(\log n/n)\) difference on these compact intervals, so (9) also holds for \(\widetilde Q_n\).

All central/calibration response queries belong to a fixed compact subset of \([-2,7]\). Take
\[
K_N=\left\lceil\frac{2\log N}{|\log\lambda|}\right\rceil+2 .
\tag{10}
\]
The population truncation error in (4) is \(O(N^{-2})\). There are \(O(K_N)\) quantile terms, so their summed nonlinear remainder is \(O_P(N^{-3/4}\log^2N)=o_P(N^{-1/2})\). The linear terms are exactly the empirical averages of the truncated quantile representers in the local proof. Their difference from the infinite representer has \(L^2\) norm \(O(\lambda^{K_N/2})=O(N^{-1})\). Thus fixed-argument profile reconstruction has a root-\(N\) expansion; this conclusion uses both a growing-series remainder bound and the influence convergence.

For moving arguments \(|x-y|\le d_N\), (8) and \(|C^kx-C^ky|\le\lambda^k|x-y|\) imply for the centered linear reconstruction
\[
O_P\!\left(N^{-1/2}\sqrt{d_N\log N}+K_N\log N/N\right).
\tag{11}
\]
Smooth coefficient factors in the quantile kernels add only \(O_P(d_N/\sqrt N)\). Sum the square-root interval terms geometrically, rather than bounding every one by the same moving-argument error. The nonvanishing logarithmic remainder is summed only \(K_N\) times. Taking \(d_N=B\log N/\sqrt N\) makes (11) \(o_P(N^{-1/2})\); (9)'s summed remainder is also negligible. This supplies the needed stochastic equicontinuity.

## B. Both response tails, the jump at 4, and interpolation

For \(U_c=A_c-W_c\) or \(W_c\), and also its first two coefficient derivatives,
\[
|U_c(z)|+|U'_c(z)|\le C(1+|z|)^m e^{-z^2/2+C|z|},
\quad z\ne4,\qquad
F'_a(z)\le C(1+|z|)^m e^{.1z^2+C|z|}.
\tag{12}
\]
These follow from the explicit Gaussian/warped-Gaussian kernels, bounded shifts, and outward lattice sums. Constants are uniform on the public coefficient rectangles and (1). The finite jump of \(W_c\) at 4 is handled separately below.

For the **step quantile**, set \(w_R(p)=U_c(\Phi^{-1}p)/\phi(\Phi^{-1}p)\) on \(|\Phi^{-1}p|\le R\), zero otherwise, and \(V_R(p)=\int_0^p w_R(t)\,dt\). Layer-cake integration on positive responses yields the exact signed-weight identity
\[
\int_{-R}^R U_c(z)[Q_n(\Phi z)-F(z)]dz
=-\int_{\mathbb R}F'(z)
 [V_R(H_n(\Phi z))-V_R(\Phi z)]dz.
\tag{13}
\]
The common total weight cancels the endpoint term. This identity is not applied directly to the polygonal quantile.

Set \(R_N=\sqrt{\tfrac32\log N}\). On event (7), the empirical and population probabilities are relatively close at the cutoff, since \(n\Phi(-R_N)=N^{1/4+o(1)}\). Segments that intersect the support of \(w_R\) correspond, eventually, to \(|z|\le R+1\); their probabilities and Gaussian densities are comparable. On smooth pieces,
\[
|w'_R(p)|=|(U'_c(z)+zU_c(z))/\phi(z)^2|
\le C(1+|z|)^m e^{z^2/2+C|z|}.
\]
Taylor expansion of (13), using (7) and \(p(1-p)\le C\phi(z)/(1+|z|)\), bounds its integrated nonlinear remainder by
\[
\frac{C\log N}{N}\int_{|z|\le R+1}
 F'(z)(1+|z|)^m e^{C|z|}dz
=N^{-.85+o(1)}.
\tag{14}
\]
There are jumps at \(-R,4,R\). For a jump \(J_0\) of \(w_R\) at \(p_0=\Phi(z_0)\), a crossing is confined to a strip of width \(O(d_0)\), where \(d_0^2=Cp_0(1-p_0)\log N/N\). Its additional remainder is bounded by
\[
C F'(z_0)|J_0|d_0^2/\phi(z_0)
\le (\log N/N)(1+|z_0|)^m e^{.1z_0^2+C|z_0|}.
\tag{15}
\]
Thus the two cutoff jumps have order (14), and the fixed jump at 4 is \(O(\log N/N)\). Omitting this kink contribution would be an error.

For the polygon replacement, on (7) and the maximum-spacing event, every adjacent uniform spacing relevant to \(|z|\le R\) is \(O(\log N/N)\), and its surrounding probability is comparable to \(\Phi(z)\) and its complement. Hence
\[
|\widetilde Q_n(\Phi z)-Q_n(\Phi z)|
\le C(\log N/N)F'(z)/\phi(z).
\]
Multiplication by \(|U_c|\) and integration gives \(N^{-.85+o(1)}\). The compact replacement error after \(K_N\) terms is \(O_P(\log^2N/N)\). This proves that the continuous interpolation does not change the first-order influence.

Population tail bias in (6) is
\[
\int_{|z|>R}|U_c(z)F(z)|dz
\le C\int_{|z|>R}(1+|z|)^m e^{-.4z^2+C|z|}dz
=N^{-.6+o(1)}.
\tag{16}
\]
The omitted linear influence has norm \(N^{-.225+o(1)}\), by the individual-source Bochner envelope \(e^{-.15z^2+C|z|}\) from the local proof; its empirical average is \(N^{-.725+o_P(1)}\). Thus removal of the cutoff is valid at the root-\(N\) scale, without interchanging two unrelated limits.

Here \(N^{o(1)}\) abbreviates a fixed bound \(C(\log N)^m e^{C\sqrt{\log N}}\); constants come from the declared envelopes and coefficient rectangles, not an unknown estimated radius or a slow diagonal selection. The compatible inequalities for \(R^2=c\log N\) are \(c<2\), \(.4c>.5\), and \(.1c<.5\). The explicit \(c=1.5\) satisfies all three.

## C. All four coefficients: selection and expansion

For pair \(i\), define the population calibration map
\[
B_i(\alpha,\beta)=
(g_1(-1+\alpha)+g_2(-h_1+\beta),\
 g_1(1+\alpha)+g_2(h_1+\beta)),\quad h_1=h(1).
\]
Its observed right-hand side is \((F_i(-1),F_i(1))\). On the public rectangle of radius \(r=.02\), the difference between any two calibration values equals an averaged derivative matrix times their coefficient difference. Its columns can be averaged separately; they obey
\[
\det\overline M_i\ge D_i:=
e^{\alpha_{i0}+\beta_{i0}-2r}
[(1-.01)^2e^{h_1-1}-(1+.01)^2e^{1-h_1}]>0.
\tag{17}
\]
To justify the averaged determinant, factor the common positive quantities $\overline A=\int_0^1e^{\alpha_t}dt$ and $\overline B=\int_0^1e^{\beta_t}dt$ from its two columns. Both rows of each column share its averaging path. The remaining positive and negative products have factors $(.99)^2$ and $(1.01)^2$; $\overline A\overline B\ge e^{\alpha_{i0}+\beta_{i0}-2r}$ gives (17). Independent unrelated entry bounds would be too crude for this constant. Since \(h_1-1=\sqrt2/10\), the bracket is strictly positive. Every entry is at most
\(B_*=(1.01)e^{h_1+1}\), so
\(\|B_i(c)-B_i(c')\|_2\ge \gamma\|c-c'\|_2\),
where the public constant \(\gamma=\min_iD_i/(2B_*)>0\).
This proves uniqueness and quantitative separation throughout the search rectangle, not merely an invertible Jacobian at the truth.

Replace profiles and right-hand sides by their reconstructed data versions. On a deterministic rectangular grid of mesh at most \(N^{-1}\), select the first lexicographic minimizer of the Euclidean residual. A numerical objective evaluated uniformly to error at most \(N^{-2}\) gives an equally valid approximate selection. There is no choice of a root nearest the unknown truth.

The crude uniform reconstruction error is \(O_P(K_N/\sqrt N)=O_P(\log N/\sqrt N)\). Comparing the minimizer with a grid point near the truth and using (17) gives the same preliminary coefficient error. The truth is at least .01 from each boundary of the search rectangle, so this preliminary estimate lies in its interior with probability tending to one.

At the truth let \(\widehat R_i(c_i)\) denote empirical calibration residual. Part A shows it is root-\(N\) and gives its joint expansion. Use (11) on the preliminary neighborhood and Taylor expand the population map:
\[
\widehat R_i(\widehat c_i)=
\widehat R_i(c_i)-M_{i,\eta}(\widehat c_i-c_i)
+o_P(N^{-1/2}).
\tag{18}
\]
Its population quadratic remainder is \(O_P(\log^2N/N)\). To prove that the actual minimum residual is small, consider only in the proof
\(c_i^*=c_i+M_{i,\eta}^{-1}\widehat R_i(c_i)\).
Equation (18) and its nearest grid point give residual \(o_P(N^{-1/2})\). Therefore the actual measurable minimizer also has that residual. Consequently
\[
\widehat c_i-c_i=M_{i,\eta}^{-1}\widehat R_i(c_i)+o_P(N^{-1/2}).
\tag{19}
\]
These are exactly the four coefficient duals \(\chi_{\eta,j}\) in the local proof, jointly. The matrices, true coefficients, and duals are proof objects, not supplied to the estimator.

## D. Target composition and actual influence

Parts A–B hold uniformly for \(c\) in the public rectangles and the first two coefficient derivatives of the target weights. For example the absolute coefficient increment of the tail linear empirical functional is bounded by
\[
|\widehat c-c|\int \sup_\theta|\partial_c U_\theta(z)|
\frac{F'(z)}{\phi(z)}|\Delta_n(\Phi z)|\,dz.
\]
The integral is \(O_P(N^{-1/2})\): its expectation is bounded by the finite individual-source Bochner envelope. The compact integral has the same property by geometric summability. Thus substituting root-\(N\) coefficients changes the empirical fluctuation by \(O_P(N^{-1})\). Taylor expansion of the population target adds \(b_{\eta,c}^T(\widehat c-c)\) with quadratic \(O_P(N^{-1})\) error; the second derivatives are bounded by the original profile envelopes.

Combining (19) with the anchor target influence gives
\[
\psi_\eta=\psi_{A,\eta}+\sum_{j=1}^4b_{\eta,c,j}\chi_{\eta,j}.
\tag{20}
\]
This is the influence of this actual statistic. Each quantile kernel has the allocation factor \(1/\pi_a=5\); the empirical linear term is therefore exactly the \(1/N\) sum in (2), not a sum with an erroneous \(1/n\) normalization. Repeated contributions from source 0, including those in calibration, are combined in the same \(\psi_{\eta,0}\) before taking its square. No sample splitting or fictitious independence is used.

## E. Lattices, numerical approximation, and fallback

Take \(M_N=\lceil4R_N+30\rceil+2\), summing \(P\) over \(-M_N,\ldots,M_N\), and \(W\) over \(1,\ldots,M_N\) or \(0,\ldots,M_N-1\). Omitted weights on the central cell and retained tail interval are bounded by a polynomial times
\(\exp[-(M_N-5)^2/2+CM_N]=N^{-12+o(1)}\).
On the observable event \(\max Y\le N^2\), all polygonal responses are bounded by \(N^2\). Thus even a deterministic bound for the lattice error is \(N^{-10+o(1)}\). The finite periodization is used with profile **differences**, so its nonzero truncated integral is not silently treated as zero.

The statistic returns zero for fewer than two observations per group or if any response exceeds \(N^2\), and clips its final target to \([-1000,1000]\). The following is a sufficient proof-onset diagnostic, not an operational fallback:
\[
\frac N5\min\{\ell(7),\ell(R_N+1)\}\ge64(\log N)^2,\qquad
\ell(x)=\phi(x)x/(1+x^2)\le\Phi(-x).
\tag{21}
\]
It holds eventually since its second effective rank is \(N^{1/4+o(1)}\). For smaller \(N\), the polygonal statistic remains defined and is not forced to zero by this diagnostic. No useful finite-sample accuracy follows before or merely from this condition. The low-count fallback is deterministic and disappears once \(n\ge2\).

Uniformly in the original class,
\[
EY^4\le M_4:=128e^4(e^8+\sqrt5e^{40.2}),\qquad
P(\max Y>N^2)\le M_4N^{-7}.
\tag{22}
\]
Indeed \(g_j\le2e^x\), shifts are at most 1, and \(h(z)\le z+.1z^2+.05\). The Gaussian quadratic moment then gives (22). Positive finite observations occur almost surely. The target envelope in the original model is strictly below 1000: \(m_0<2\), \(m_h\le e^{.675}/\sqrt{.8}<3\), and \(2(e+1)^2(m_0+m_h)<160\). Thus clipping is inactive with probability tending to one. Bounded fallback changes (2) on an event of vanishing probability.

For numerical realization, round data and evaluate known elementary functions to absolute error \(N^{-40}\); use bisection for \(h^{-1}\) with bracket width \(N^{-40}\); use midpoint mesh at most \(N^{-6}\), splitting at \(z=4\). On the response-cap event, the polygon is globally \(N^3\)-Lipschitz in probability. Finite target integrands are piecewise \(N^{3+o(1)}\)-Lipschitz, including the compact contraction integral. Hence quadrature error is \(N^{-3+o(1)}\). Repeated contraction-query errors accumulate at most a constant times \(N^{-40}/(1-\lambda)\). All primitive, input and summation errors remain \(o(N^{-1/2})\), even with \(N^{6+o(1)}\) operations. Objective errors can be made below \(N^{-2}\) uniformly over the finite coefficient grid.

The arithmetic model is rational interval evaluation of known elementary functions, with certified enclosures refined until these absolute tolerances hold. Square roots and the monotone inverse admit rational bisection; exponential and Gaussian integrals admit Taylor series with remainder bounds after range reduction; \(\pi\) is computable to any specified rational enclosure. Polygon interpolation is continuous at knots, so interval enclosure does not require deciding an exact transcendental rank comparison. Near coefficient ties, select using objective approximations with the displayed error, rather than demanding an exact tie decision. A sufficient working-bit schedule is \(\lceil10000+200\log_2N\rceil\), with enclosure refinement if needed; this is an explicit accuracy obligation, not a claim that ordinary floating-point arithmetic certifies it. See implementation limitations in [estimator.md](estimator.md).

## Explicit simultaneous schedule and scope

| Item | Fixed choice | Contribution at root-\(N\) comparison |
|---|---|---|
| Contraction | (10) | population \(O(N^{-2})\); summed Bahadur \(O_P(N^{-3/4}\log^2N)\) |
| Cutoff | \(R_N^2=1.5\log N\) | population \(N^{-.6+o(1)}\); linear tail \(N^{-.725+o_P(1)}\) |
| Tail nonlinear/interpolation error | specified polygon, no fitted smoothing | \(N^{-.85+o_P(1)}\) |
| Central interpolation | knot spacing \(1/n\) | \(O_P(\log^2N/N)\) |
| Coefficient grid | mesh \(\le N^{-1}\) in public rectangles | \(o_P(N^{-1/2})\) residual; (19) |
| Lattice | \(M_N=\lceil4R_N+30\rceil+2\) | \(N^{-10+o(1)}\) under cap |
| Quadrature | mesh \(\le N^{-6}\), split at 4 | \(N^{-3+o(1)}\) |
| Inverse/primitive/input error | \(N^{-40}\) | negligible after all operations |
| Residual objective accuracy | \(\le N^{-2}\) | negligible for (19) |
| Response-cap fallback | \(Y>N^2\) | probability \(\le M_4N^{-7}\) |

The previously proved fixed-path DQM supplies LAN for the balanced product experiment. Apply the joint CLT to the influence sum and the fixed-path score sum, then the likelihood change-of-measure lemma. Their covariance is \(DI_\eta[v]\) by the full score equation; the mean shift cancels the local change in \(I\). Contiguity transfers the \(o_P\) remainder. This proves regularity along the specified \(N^{-1/2}\) paths, without claiming a uniform nonlinear expansion over unrestricted moving directions.

**OPEN / not claimed:** efficient attainment, a consistent implementable variance estimator, studentized coverage, \(O(1/N)\) MSE, a neighborhood-uniform CLT, a global Hellinger modulus, and minimax optimality. A CLT and a bounded final statistic alone do not give the needed uniform integrability for an MSE rate. A variance estimator based on (20) would need justified estimates of its unknown profile/density derivatives or a separate derivative-free construction; neither is presumed here.

The theorem advances the saved result from oracle influence existence to a specified, source-only, pointwise regular estimator. It does not establish a new general learning principle, practical superiority, or biological validity.
