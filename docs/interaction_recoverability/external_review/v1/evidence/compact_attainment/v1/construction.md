# Compact construction, same five-law model

**PROVED HERE.** The compact anchor equation and complete finite statistical construction below retain the saved model. Two positive analytic profiles satisfy \(g_j(0)=1\), \(g_j(-\infty)=0\). At each fixed truth,
\[
\max_{j,k\le3}\sup_x e^{-x}|g_j^{(k)}(x)-e^x|\le .01.
\]
The original derivative envelopes, admissible analytic paths and known \(Z\sim N(0,1)\) remain as in the [original model](../../unknown_profiles/model.md) and [saved theorem](../../local_attainment/v1/theorem.md). The known warp is
\(h(z)=z+.1z\sqrt{1+z^2}\).
The unknown coefficient centers are publicly specified as
\(\alpha_0=(2/3,3/4)\), \(\beta_0=(3/5,4/5)\); the truth is within .01 of each center and the public search rectangles have radius .02.

The observations are \(n\) independent responses from each of \(0,e_1,e_2,e_3,e_4\), \(N=5n\). Their increasing response maps are
\[
F_0(z)=g_1(z)+g_2(hz),\quad
F_i(z)=g_1(z+\alpha_i)+g_2(hz+\beta_i),\ i=1,2,
\]
\[
F_3(z)=g_1(z+1)+g_2(hz),\quad F_4(z)=g_1(z)+g_2(hz+1).
\]
The unit private loadings are known. Neither latent realizations, response derivatives, profiles, target labels nor true coefficients are observed. Quantile queries at transformed arguments evaluate the same five source laws; they add no observed actions.

## Population and tangent algebra

Put \(T_jx=h^{-1}(hx+j)\). For any positive integers \(p,q\),
\[
E_{pq}(x)=F_0(T_qx)-F_0(x)
-\sum_{j=0}^{q-1}(F_4-F_0)(T_jx)
-\sum_{k=1}^{p}(F_3-F_0)(T_qx-k).
\]
The first sum is \(g_2(hx+q)-g_2(hx)\). Subtracting it from the first difference leaves \(g_1(T_qx)-g_1(x)\). The second sum is \(g_1(T_qx)-g_1(T_qx-p)\). Thus
\[
E_{pq}(x)=g_1(T_qx-p)-g_1(x). \tag{1}
\]
Only \(p=4,q=5\) is used henceforth.

For an admissible profile path with derivative \(r=(r_1,r_2)\), the anchor tangent curves are
\(q_0=r_1+r_2\circ h\),
\(q_3=r_1(\cdot+1)+r_2\circ h\),
\(q_4=r_1+r_2(h+1)\).
Repeating the two finite telescopes gives
\(E_{45}(q;x)=r_1(Cx)-r_1(x)\).
The four coefficient directions have zero derivatives in these three anchor laws, so the same identity holds along mixed paths. This is an actual tangent calculation, not differentiation of an empirical fitted profile.

Let \(C=T_5-4\), \(J=[-7/4,-3/4]\), \(b=-7/4\), \(r=-1\), \(\lambda=271/290\). The endpoint certificates are
\[
h(9/4)-h(-7/4)=4+(9\sqrt{97}+7\sqrt{65})/160<4+153/160<5,
\]
\[
h(13/4)-h(-3/4)=4+(13\sqrt{185}+15)/160>4+184/160>5.
\]
Therefore \(T_5(J)\subset(9/4,13/4)\), and \(C(J)\subset J\).
Since \(h'\) is positive, even and increasing in \(|x|\),
\[
h'(7/4)=1+57/(20\sqrt{65})<271/200,\quad
h'(9/4)=1+89/(20\sqrt{97})>29/20.
\]
Squaring positive quantities reduces the strict bounds to
\(570^2<71^2\,65\) and \(89^2>81\,97\).
Consequently \(0<C'<\lambda\). Every \(T_jx\) lies between \(x\) and \(T_5x\), and every \(T_5x-k\) lies between \(Cx\) and \(T_5x-1\). All source queries are in
\[
D_{\rm compact}=[-7/4,13/4]. \tag{2}
\]
After cancelling \(F_0(x)\), \(E_{45}\) has 18 signed response evaluations (nine positive, nine negative), compared with six for the old anchor equation.

## Normalized reconstruction and calibration charts

Define
\[
D_K(F;x,y)=-\sum_{k=0}^{K-1}[E_{45}(F;C^kx)-E_{45}(F;C^ky)].
\]
Exactly,
\[
D_K=g_1(x)-g_1(y)-\{g_1(C^Kx)-g_1(C^Ky)\}. \tag{3}
\]
The last term is at most \(1.01e^{-3/4}\lambda^K|x-y|\) on \(J\).
Normalization is
\[
g_1(r)=1-(F_3-F_0)(-1). \tag{4}
\]
For arbitrary \(x\), choose the deterministic integer \(m(x)=\lceil b-x\rceil\), \(t=x+m\in[b,b+1]\) (with the left endpoint convention at integers). Put
\[
\widehat g_{1,K}(x)=1-(\widetilde F_3-\widetilde F_0)(r)+D_K(\widetilde F;t,r)+H_m(\widetilde F;x),
\]
where \(H_m=-\sum_{j=0}^{m-1}(\widetilde F_3-\widetilde F_0)(x+j)\) for \(m>0\);
\(H_m=\sum_{k=1}^{-m}(\widetilde F_3-\widetilde F_0)(x-k)\) for \(m<0\); and \(H_0=0\).
Set \(\widehat g_{2,K}(y)=\widetilde F_0(h^{-1}y)-\widehat g_{1,K}(h^{-1}y)\).
These are finite response operators, without estimating any derivative.

For calibration, \(B_i(c)\) has rows
\(g_1(-1+\alpha)+g_2(-h_1+\beta)\) and
\(g_1(1+\alpha)+g_2(h_1+\beta)\), \(h_1=h(1)\).
Throughout both public rectangles:

| Query | Original interval enclosing query | Integer \(m\) | Transported interval |
|---|---|---:|---|
| \(-1+\alpha\) | \((-.36,-.23)\) | -1 | \((-1.36,-1.23)\) |
| \(1+\alpha\) | \((1.64,1.77)\) | -3 | \((-1.36,-1.23)\) |
| \(h^{-1}(-h_1+\beta)\) | \((-.6,-.25)\) | -1 | \((-1.6,-1.25)\) |
| \(h^{-1}(h_1+\beta)\) | \((1.4,1.7)\) | -3 | \((-1.6,-1.3)\) |

Here \(1.14<h_1<1.15\), \(.58\le\beta\le.82\), and evaluating the monotone \(h\) at the displayed rational endpoints proves the inverse bounds. Every transported query is at least .15 from a join. The additional finite translations and normalizer remain inside (2). Moving-query equicontinuity is invoked only on these fixed charts.

## Signed target and new outward identity

Write \(\rho_1=\phi\), \(\rho_2(y)=\phi(h^{-1}y)/h'(h^{-1}y)\), and
\[
K_{j,c}(x)=\rho_j(x-c_{j1}-c_{j2})-\rho_j(x-c_{j1})-\rho_j(x-c_{j2})+\rho_j(x).
\]
Then the biological-style joint-action contrast in this mathematical model is
\[
I(g,c)=\sum_{j=1}^2\int K_{j,c}(x)g_j(x)\,dx.
\]
With \(A_c(z)=h'(z)K_{2,c}(h(z))\) and \(L_c=K_{1,c}-A_c\),
\[
I=\int A_c F_0+\int L_c g_1,\quad
P_c(t)=\sum_{j\in\mathbb Z}L_c(t+j),
\]
\[
W_{b,c}(z)=
\begin{cases}
\sum_{j\ge1}L_c(z+j),&z\ge b,\\
-\sum_{j\ge0}L_c(z-j),&z<b.
\end{cases}
\]
Partition the real line into \(t+j\), \(t\in[b,b+1]\), expand \(g_1(t+j)\) using forward/backward unit differences, and sum their absolutely integrable coefficients. The result is
\[
I=\int_{\mathbb R}[(A_c-W_{b,c})F_0+W_{b,c}F_3]\,dz
+\int_b^{b+1}P_c(t)[g_1(t)-g_1(r)]\,dt. \tag{5}
\]
The sign of the \(F_3\) term is positive. The jump is
\(W(b+)-W(b-)=P_c(b)\), at \(b=-7/4\), not at the old point 4.
All sums/integrals converge absolutely by the Gaussian envelopes in the theorem. The infinite periodization has zero integral because \(\int L_c=0\). Its finite version generally does not.

Define \(P_{M,c}=\sum_{j=-M}^M L_c(t+j)\); truncate the two \(W\) sums to \(1,\ldots,M\) and \(0,\ldots,M-1\), respectively. The finite ideal target functional is
\[
H_{KMR}(\widetilde F,c)=
\int_{-R}^{R}[(A_c-W_{M,c})\widetilde F_0+W_{M,c}\widetilde F_3]\,dz
+\int_b^{b+1}P_{M,c}(t)D_K(\widetilde F;t,r)\,dt. \tag{6}
\]
In particular the finite integral uses the difference \(D_K\). Replacing it by the reconstructed profile and dropping its constant would add the nonzero quantity
\(\widehat g_1(r)\int_b^{b+1}P_{M,c}\), an error.

At population curves, the exact difference between (5) and (6) consists of: the omitted response-tail integral; \(\int_{-R}^R(W-W_M)(F_3-F_0)\); \(\int_J(P-P_M)(g_1-g_1(r))\); and
\(\int_J P_M(t)[g_1(C^Kt)-g_1(C^Kr)]dt\).
These explicit remainders also fix every sign.

## Comparable central-operator upper bounds

Use the stratified Hilbert norm \(\|v\|^2=\sum_a(1/5)E v_a^2\). For \(D=[l,u]\), \(B=\max(|l|,|u|)\), define public envelopes
\[
F_1=1.01[e^{u+1}+h'(B)e^{h(u)+1}],\quad
F_2=1.01[e^{u+1}+(.2+h'(B)^2)e^{h(u)+1}],
\]
\[
C_D=\sqrt5\left\{\frac{\sqrt{u-l}(F_2+BF_1)}{2\phi(B)}
+\frac{F_1\sqrt{\phi(0)}}{\phi(B)}\right\}.
\]
The individual quantile representers in the theorem obey
\(\|\kappa_{a,z}-\kappa_{a,w}\|\le C_D\sqrt{|z-w|}\).
Indeed their centered-indicator increment has norm at most
\(\sqrt{\phi(0)|z-w|}\); their smooth coefficient derivative is bounded by
\((F_2+BF_1)/\phi(B)\). Use \(|z-w|\le\sqrt{u-l}\sqrt{|z-w|}\).

For an equation with \(m\) signed terms and query maps Lipschitz at most \(L\), the infinite central difference has bound
\[
\|d_{x,y}\|\le
\mathcal B_D\sqrt{|x-y|},\qquad
\mathcal B_D=\frac{m C_D\sqrt L}{1-\sqrt\lambda}. \tag{7}
\]
The old construction admits \(D=[3,7],m=6,L=2.3,\lambda=575/613\).
The new admits \(D=[-7/4,13/4],m=18,L=271/200,\lambda=271/290\).
This compares the central difference operators; normalization and calibration are separate terms.

**NUMERICALLY CHECKED:** these same conservative formulas give old
\(7.3152316001\,10^{20}\) and new \(1.7921492674\,10^9\).
Their ratio is \(2.4498872563\,10^{-12}\). Neither number is the actual operator norm, an efficient variance, or a lower bound. Comparing two loose upper bounds does not establish a variance improvement. Both infinite operators represent the same on-model target derivative, but can differ on the orthogonal complement of the model tangent space.

