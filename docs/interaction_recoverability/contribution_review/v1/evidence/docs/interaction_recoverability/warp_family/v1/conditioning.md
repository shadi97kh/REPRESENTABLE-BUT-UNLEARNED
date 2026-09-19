# Explicit conditioning and a complete-experiment lower bound

All statements retain [problem.md](problem.md): known \(\theta\in[1/2,2]\), \(0<\delta\le1/3{,}200{,}000\), full original profile envelopes, original coefficient box, and \(N=5n\). All bounds below are derived analytically. No numerical values were generated.

## 1. Warp and moment envelopes

**PROVED HERE.** Write
\[
\omega=\delta\sqrt\theta,\qquad a=\delta/(2\sqrt\theta),\qquad
\pi_*=\min_a\pi_a=1/5.
\]
Here the un-subscripted \(\pi\) inside Gaussian constants is the usual circle constant. The inequality
\[
|x|\sqrt{1+\theta x^2}\le\sqrt\theta x^2+\frac1{2\sqrt\theta}
\]
follows by squaring positive sides, or by \(\sqrt{t^2+1}\le t+(2t)^{-1}\). Thus
\[
|h(z)|\le|z|+\omega z^2+a,\quad
h'(z)\le1+\delta+2\omega|z|,\quad |h''(z)|\le2\omega. \tag{1}
\]
For \(p>0\) with \(2p\omega<1\), a public bound is
\[
\overline m_{h,p}
=\frac{2e^{pa}}{\sqrt{1-2p\omega}}
\exp\!\left\{\frac{p^2}{2(1-2p\omega)}\right\}
\ge E e^{p h(Z)}. \tag{2}
\]
Bound \(e^{p|Z|}\) by \(e^{pZ}+e^{-pZ}\) and integrate each Gaussian quadratic exponential to obtain (2). The factor 2 is conservative.
For actual sources, \(g_j\le2e^x\) and shifts are at most 1. Therefore
\[
E Y_a^4\le M_4:=128e^4(e^8+\overline m_{h,4}),\qquad
P(\max_{a,i}Y_{ai}>N^2)\le M_4N^{-7}. \tag{3}
\]
The verified interval has \(\omega<1/100\), so the required fourth-moment condition \(8\omega<1\) holds with a fixed margin. These exponents are derived for this family; the old \(\delta=.1,\theta=1\) exponents are not copied.

## 2. Quantile kernels and compact inversion

**PROVED HERE.** At any fixed truth, the centered source-specific representer is
\[
\kappa_{a,z}(y)=\frac{F_a'(z)}{\pi_a\phi(z)}
\{\Phi(z)-\mathbf1[y\le F_a(z)]\},
\]
in coordinate \(a\), zero in the other coordinates. Its squared stratified norm is
\[
\|\kappa_{a,z}\|_\pi^2
=\frac{F_a'(z)^2\Phi(z)(1-\Phi(z))}{\pi_a\phi(z)^2}.
\]
Using (1), the slightly enlarged shift bound \(1.1\), and
\(\Phi(z)(1-\Phi(z))\le\frac12e^{-z^2/2}\), gives
\[
\|\kappa_{a,z}\|_\pi\le Q(z):=
2e^{1.1+a}\sqrt{\frac{\pi}{\pi_*}}
(2+\delta+2\omega|z|)
e^{(1/4+\omega)z^2+|z|}. \tag{4}
\]
This bound retains Gaussian extreme-quantile amplification.

On the central query interval \([-2,2]\), \(|h(z)|\le3\), \(h'\le2\), \(|h''|\le1\). Indeed \(|s(2)|\le2M\) and \(2\delta M<1\) on the proved interval. Source derivatives obey
\[
F'\le F_1:=6e^4,\qquad |F''|\le F_2:=52e^4.
\]
The original bound on \(g''\), not an exponential-profile identity, supplies \(F_2\). Define
\[
C_\kappa=
\frac{F_2+2F_1+F_1\sqrt{\phi(0)}}{\sqrt{\pi_*}\phi(2)}.
\]
Decomposing a kernel difference into its smooth coefficient change and a centered indicator increment gives
\[
\|\kappa_{a,x}-\kappa_{a,y}\|_\pi\le C_\kappa\sqrt{|x-y|}.
\]
For the smooth term use \(|(F'/\phi)'|\le(F_2+2F_1)/\phi(2)\), indicator standard deviation at most \(1/2\), and \(|x-y|\le2\sqrt{|x-y|}\). For the indicator increment use variance at most \(|\Phi(x)-\Phi(y)|\le\phi(0)|x-y|\).

Each of the six commutator maps is 2-Lipschitz. With
\[
\lambda=1-\kappa\delta/2,\quad \kappa=1/64,\quad
C_E=6\sqrt2 C_\kappa,\quad
D_c=\frac{C_E}{1-\sqrt\lambda}\le\frac{4C_E}{\kappa\delta}, \tag{5}
\]
the infinite central difference influence \(d_{x,r}\) has norm at most \(D_c\sqrt{|x-r|}\). For \(r=-1,x\in[-2,0]\), adding the normalizer gives the profile-evaluation bound
\[
G_c=2Q(2)+D_c. \tag{6}
\]
The series tail after \(K\) terms has norm at most \(D_c\lambda^{K/2}\). This is convergence in the direct sum of individual source spaces, not just convergence of population values. In particular an ordinary source \(L^2\) norm is at most the corresponding stratified bound divided by \(\sqrt{\pi_a}\).

These statements bound central profile evaluations and, by integration, compact signed profile functionals. They are not the complete interaction variance or an efficient bound.

## 3. Transports and the complete calibration cost

**PROVED HERE.** Use the radius \(R\) and fixed charts from [commutator_and_calibration.md](commutator_and_calibration.md). Put
\[
D=R+3,\quad L=\lceil R+4\rceil,\quad
\gamma=\tfrac12e^{.4}-2e^{-1.9}.
\]
Let \(U_x\) denote the profile-1 influence with its prescribed fixed chart, and
\(V_y=\kappa_{0,h^{-1}y}-U_{h^{-1}y}\) the profile-2 influence. At the calibration arguments,
\[
\|U_x\|_\pi\le G_c+2LQ(D).
\]
For \(\sigma=-,+\), \(z_\sigma=\sigma R\), \(w_-=e^R,w_+=e^{-h(R)}\), define the scaled calibration residual representers
\[
\zeta_{i,\sigma}
=w_\sigma\{\kappa_{i,z_\sigma}
-U_{z_\sigma+\alpha_i}
-V_{h(z_\sigma)+\beta_i}\}.
\]
Both row norms are bounded by
\[
Z_R=e^R[\,2G_c+(4L+2)Q(D)\,]. \tag{7}
\]
The bound uses \(w_\sigma\le e^R\); it is deliberately conservative for the positive row. The actual scaled Jacobian \(J_i\) has inverse infinity norm at most \(\gamma^{-1}\). Therefore each component of
\(\chi_i=J_i^{-1}(\zeta_{i,-},\zeta_{i,+})^\top\) satisfies
\[
\|\chi_j\|_\pi\le Z_R/\gamma. \tag{8}
\]
The algebraic separation margin stays positive over the full profile class. The statistical cost in \(Q(D)\) and \(e^R\) does not disappear.

## 4. An elementary complete-influence upper bound

**PROVED HERE.** Use the unit cell \([-1,0]\), the signed target kernels \(K_1,K_2\), \(A=h'K_2\circ h\), \(\mathcal L=K_1-A\), and outward weight \(W\), as specified in [attainment.md](attainment.md). A four-term signed sum of densities has \(L^1\) norm at most 4. Consequently
\[
\int|\mathcal L|\le8,\qquad
\int_{-1}^0|P(t)|\,dt\le8.
\]
The central target influence has norm at most \(8D_c\).

For an explicit outer bound, put
\[
B_0=11/5,\quad S=11e^{(B_0+1)^2},\quad
\nu=\tfrac14-\omega,\quad A_*=2+\delta,\quad C_*=2\omega.
\]
If \(w=h^{-1}(h(z)-b)\), \(0\le b\le B_0\), then \(|w-z|\le b\), and \(\phi(w)\le\phi(0)e^{-z^2/2+B_0|z|}\). The bounds \(h'\ge1\) and (1) give, separately for \(A\) and \(\mathcal L\),
\[
|A(z)|,\ |\mathcal L(z)|
\le4\phi(0)(A_*+C_*|z|)e^{-z^2/2+B_0|z|}.
\]
For the outward sums from \(-1\), the same envelope bounds \(W\) after multiplication by \(S\). To verify a public scalar bound for the sum, use
\[
\sum_{j\ge0}(1+j)e^{-j^2/2+(B_0+1)j}
\le e^{(B_0+1)^2}\sum_{j\ge0}(1+j)e^{-j^2/4}
<11e^{(B_0+1)^2}.
\]
Here \(\sum e^{-j^2/4}\le1+\sqrt\pi\), and \(j\le2e^{j^2/8}\) gives
\(\sum j e^{-j^2/4}\le2(1+\sqrt{2\pi})\). The extra 1 in \(B_0+1\) covers \(z\in[-1,0]\) on the positive branch. Polynomial prefactors are bounded using \(C_*\le A_*\).

Multiplying by (4), then using
\((B_0+1)|z|\le\nu z^2/2+(B_0+1)^2/(2\nu)\), proves
\[
\begin{split}
D_t:={}&8e^{1.1+a}\phi(0)\sqrt{\frac{\pi}{\pi_*}}(1+2S)
e^{(B_0+1)^2/(2\nu)}\\
&\times\left[
A_*^2\sqrt{\frac{2\pi}{\nu}}
+\frac{4A_*C_*}{\nu}
+\frac{C_*^2\sqrt{2\pi}}{\nu^{3/2}}
\right]
\end{split} \tag{9}
\]
is an upper bound for the outer influence norm. These are the elementary integrals of \((A_*+C_*|z|)^2e^{-\nu z^2/2}\). The sum of absolute source-vector norms is integrable before any source covariance is combined.

At true coefficients in \([.5,1]^4\), the complete target coefficient derivative obeys
\[
\|b_c\|_1\le B_I:=4(e^2+e)(m_0+\overline m_{h,1}),\qquad m_0=e^{1/2}.
\]
For example \(|\partial_{\alpha_1}I|\le2(e^2+e)m_0\), using \(g_1'\le2e^x\), with the analogous second-profile bound. Combining (5), (7)–(9) gives
\[
\boxed{\|\psi_\eta\|_\pi\le
\Psi(\delta,\theta):=D_t+8D_c+B_I Z_R/\gamma.}\tag{10}
\]
Every constant is public and elementary. For each source separately,
\[
E_a\psi_{\eta,a}^2\le\Psi(\delta,\theta)^2/\pi_a,
\quad V_\eta=\|\psi_\eta\|_\pi^2\le\Psi(\delta,\theta)^2.
\]
The variance is that of the constructed influence, not the minimum-norm influence.

## 5. Dependence on inverse separation

**PROVED HERE — sufficient orders, not sharpness.** From the known radius formula,
\[
R^2=\frac{3}{\delta\sqrt\theta}-\frac1{2\theta}+O(\delta),
\]
uniformly in \(\theta\in[1/2,2]\). In (10), \(D_t,B_I,\gamma^{-1}\) remain bounded over the proved interval; \(D_c=O(\delta^{-1})\) and \(L=O(\delta^{-1/2})\). The dominating conservative factor is \(e^R Q(R+3)\). Since \(\omega R^2=O(1)\),
\[
\log(1+\Psi)\le
\frac{3}{4\delta\sqrt\theta}
+O(\delta^{-1/2})+O(\log(1/\delta)), \tag{11}
\]
\[
V_\eta\le\Psi^2\le
\exp\!\left\{\frac{3}{2\delta\sqrt\theta}
+O(\delta^{-1/2})+O(\log(1/\delta))\right\}. \tag{12}
\]
The \(O\)-constants here can be chosen uniformly in \(\theta\); (10) is the explicit bound without \(O\)-notation.

This is exponential, not polynomial, sufficient deterioration for this construction. It is not an exponential lower bound. Compact central queries do remove the old moving central anchor, but they do not remove calibration queries at \(\pm R\), or their finite transports. Gaussian tail probabilities satisfy \(\Phi(-R)\) of order \(R^{-1}e^{-R^2/2}\); even reaching the ordinary tail-quantile asymptotic regime requires enormous counts under these sufficient small-\(\delta\) constants. No useful finite-sample regime is demonstrated.

The depth schedule obeys
\[
K_N=\left\lceil\frac{2\log N}{|\log(1-\kappa\delta/2)|}\right\rceil+2
\le\left\lceil\frac{4\log N}{\kappa\delta}\right\rceil+2.
\]
A putative uniform transfer along \(\delta_N\downarrow0\) cannot drop the dependence in \(K_N\), \(\phi(R+4)^{-1}\), or the inverse-quantile derivative constants. In particular \(n\Phi(-(R+4))\to\infty\) and the scaled accumulated quantile remainder would need separate control. The pointwise proof supplies neither uniformly. This is a failure of that proposed uniform inference from the proof, not an impossibility result for all estimators.

For example, with the explicit \(C_D\) in attainment Section 3, this proof would need at least the sufficient remainder condition
\[
e^R C_D(K_N+L+1)N^{-1/4}\log N\longrightarrow0
\]
to carry its scaled compact quantile error through the root-\(N\) expansion, as well as the other series/tail conditions. Known \(R\) grows with inverse separation, and this displayed inequality is not implied by \(N\delta_N^2\to\infty\). No such uniform theorem has been proved here.


## 6. A rederived lower bound for all five training laws

**PROVED HERE / STANDARD RESULT APPLIED.** The lower-bound experiment is a genuine subclass of the original profile and coefficient class. Set \(g_1=g_2=e^x\), \(m=\sqrt e\), \(D_0=e-m\), and use amplitudes
\[
(u_1,v_1)=(m,e),\qquad (u_2,v_2)=(t,m+e-t),\quad t\in[m,e].
\]
Here \(u_i=e^{\alpha_i},v_i=e^{\beta_i}\). This varies coefficients, not profile normalization. Every logarithmic coefficient remains in \([.5,1]\). Only source 2 changes; the reference, source 1 and both anchors are identical across alternatives.

For its actual density score along \(t\), let \(u=t,v=m+e-t\), \(b=1+\delta s'(z)\), \(a_z=\delta s(z)\). The response velocity divided by its derivative is
\[
f(z)=\frac{1-e^{a_z}}{u+vb e^{a_z}},\qquad S(z)=zf(z)-f'(z).
\]
Since \(u,v\ge m\), \(b\ge1\),
\[
|f|\le\frac{\delta|s|}{2m},\qquad
|\partial_{a_z}f|\le\frac1{2m},\qquad
|\partial_b f|\le\frac1m.
\]
For the first inequality use \(|1-e^a|/(1+e^a)\le |a|/2\). For the second maximize \(w(u+vb)/(u+vbw)^2\), \(w>0\), obtaining \(1/(4u)+1/(4vb)\). For the third use \(|\partial_bf|=(ve^a/(u+vbe^a))|f|\le|f|/b\le1/m\).
The last step uses the additional direct bound \(|f|\le\max\{1/u,1/(vb)\}\le1/m\), not the earlier \(\delta|s|/(2m)\) bound.
Uniformly in the warp family,
\[
|s(z)|\le2z^2+1,\quad s'(z)\le1+3|z|,\quad |s''(z)|\le3.
\]
Consequently
\[
|S(z)|\le\frac{\delta}{m}(|z|^3+2|z|+3.5)
\le\frac{7\delta}{m}(1+|z|^3),
\quad ES^2\le1568\delta^2/m^2<1600\delta^2/m^2. \tag{13}
\]
The Gaussian moment bound used here is \(E(1+|Z|^3)^2\le2(1+EZ^6)=32\). The historical numerical constant is not reused.

Let \(m_h=Ee^{h(Z)}\). Oddness and \(h(z)\ge z\) for \(z\ge0\) imply \(m_h=E\cosh(h(Z))\ge m_0\). The target derivative along the whole segment has magnitude
\[
|I'(t)|=|m_0(m-1)-m_h(e-1)|\ge k_0:=m_0D_0>0. \tag{14}
\]
Use \(H(P,Q)^2=\int(\sqrt p-\sqrt q)^2\). Integrating the square-root-density derivative along a segment of length \(\eta\), with the uniform score bound (13), gives
\[
H(P_{t+\eta,2},P_{t,2})\le20\eta\delta/m.
\]
This controls the entire nonlinear likelihood path, rather than assuming its linear approximation. Product affinity gives \(H(P^{\otimes n},Q^{\otimes n})\le\sqrt n H(P,Q)\), and total variation is at most \(H\). Choose
\[
t_-=m,\quad t_+=m+\eta,\qquad
\eta=\min\{D_0,m/(40\delta\sqrt n)\}.
\]
The **complete five-group** experiment then has TV at most \(1/2\). The ordinary two-point testing reduction for squared loss yields
\[
\boxed{\inf_{\widehat I}\sup_{\eta\ {\rm in\ original\ class}}
E_\eta(\widehat I-I_\eta)^2
\ge\frac{k_0^2}{16}
\min\left\{D_0^2,\frac{m^2}{1600n\delta^2}\right\}.}\tag{15}
\]
The supremum here is over truths, with \(\delta,\theta\) fixed and known. In terms of total \(N\), replace \(1/n\) with \(5/N\).
All five source laws are accounted for; the other four carry no distinction in these particular alternatives. No Wasserstein-to-indistinguishability implication is used.

This is a global subclass lower bound of order \(\min\{1,(N\delta^2)^{-1}\}\), not a pointwise risk bound at every truth. At an interior exponential reference, the related coefficient direction \(d_{\alpha_2}=e^{-\alpha_2},d_{\beta_2}=-e^{-\beta_2}\) has score norm \(O(\delta)\) and target derivative tending to \(m_0(e^{\alpha_1}-e^{\beta_1})\ne0\); it also gives the familiar local information lower order \(\delta^{-2}\) when these first-pair coefficients differ.

**OPEN.** The polynomial lower order and exponential constructed-influence upper order are unmatched. Moreover, (12) is a variance bound for an attained pointwise CLT, not a finite-sample MSE upper bound. The record does not establish sharp minimax dependence, necessity of extreme calibration, efficient variance, or a useful uniform shrinking-separation regime.
