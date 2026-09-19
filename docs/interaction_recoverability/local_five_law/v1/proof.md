# Local five-law regularity at the specified interior reference — v1

**Decision: local regular. PROVED HERE** means proved below for the stated model, not independently established publication novelty. The result is existence of a finite-norm influence function for every admissible analytic direction, including all four coefficients. It is not an attainment theorem for an estimator. The earlier [global unresolved decision](../../direct_target/v1/decision.md) is unchanged.

## 1. Experiment, directions, and score domain

Use the original [unknown-profile model](../../unknown_profiles/model.md), with known latent law $Z\sim N(0,1)$, known warp

\[
h(z)=z+\tfrac1{10}z\sqrt{1+z^2},\quad
\alpha=(2/3,3/4),\quad\beta=(3/5,4/5),\quad g_1=g_2=e^x.
\]

The five independent strata are $a=0,1,2,3,4$, denoting $0,e_1,e_2,e_3,e_4$, with asymptotic proportions $\pi_a=1/5$. There is no observed latent $Z$, paired outcome, or joint-action training label. Let

\[
A=(0,\alpha_1,\alpha_2,1,0),\qquad
B=(0,\beta_1,\beta_2,0,1),\quad
F_a(z)=e^{z+A_a}+e^{h(z)+B_a}.
\]

Both profiles remain unknown. Admissible directions are $v=(r_1,r_2,d)$, where each $r_j$ is real analytic, $r_j(0)=0$, and $|r_j^{(k)}(x)|\le C_v e^x$, $k=0,1,2,3$. Every $d\in\mathbb R^4$ is allowed, ordered as $\alpha_1,\alpha_2,\beta_1,\beta_2$. The path $g_{j,t}=e^x+tr_j(x)$, $c_t=c+td$ is in the original class on a nonzero two-sided interval: take $|t|C_v\le1/4$ and keep each coefficient inside its box. Then $g'_{j,t}\in[3e^x/4,5e^x/4]$, $|g''_{j,t}|,|g'''_{j,t}|\le5e^x/4<10e^x$; normalization and the negative-infinity limit persist. Integrating the positive derivative gives positivity. The interval may depend on the direction.

Use the original tangent norm

\[
\|v\|_X^2=|d|^2+
\sum_{j=1}^2\sum_{k=0}^3\int_{\mathbb R}
  |r_j^{(k)}(x)|^2\frac{e^{-2x}}{1+x^2}\,dx.
\]

Here $X$ denotes the completion of these actual analytic directions in this norm. No smooth compact-support space or finite dictionary replaces them. Let $\mathcal H=\bigoplus_a L^2_0(P_a)$, with inner product $\sum_a\pi_a E_a[p_aq_a]$, and identify each component with its pullback under $F_a(Z)$.

Write $q_a=\dot F_a$. Precisely,

\[
q_a(z)=r_1(z+A_a)+r_2(h(z)+B_a)
 +e^{z+A_a}(a_1d_{\alpha_1}+a_2d_{\alpha_2})
 +e^{h(z)+B_a}(a_1d_{\beta_1}+a_2d_{\beta_2}),
\quad u_a=q_a/F'_a,
\quad S_a(v)=zu_a-u'_a.
\tag{1}
\]

**STANDARD RESULT APPLIED, with domain checked here.** Differentiating the monotone change-of-variable density gives (1). For every original direction $u_a,u'_a$ have polynomial envelopes, $u_a\phi\to0$, and Gaussian integration by parts gives

\[
E S_a=0,\qquad \int S_a^2\phi=\int(u_a^2+(u'_a)^2)\phi.
\tag{2}
\]

For example, expand $(zu-u')^2$ and use $\int z(u^2)'\phi=\int(z^2-1)u^2\phi$. There are no discarded boundary terms. Formula (1) defines a bounded map $S:X\to\mathcal H$: write $r_j^{(k)}=e^x f_{jk}$ in the score. The remaining coefficients are polynomially bounded in $z$, after division by $F'_a\ge e^{z+A_a}+e^{h(z)+B_a}$. For the first profile, change variable $x=z+A_a$; for the second, $x=h(z)+B_a$. The resulting weights, including polynomial factors, are bounded by constants times $(1+x^2)^{-1}$: Gaussian decay for the first and exponential-in-$|x|$ decay for the second. This bounds the integrals in (2) by the $k=0,1$ terms of $X$.

The same monotone paths are differentiable in quadratic mean. To verify applicability of [P3](../../unknown_profiles/proofs.md), differentiate $\sqrt{p_t(y)}$ as $f_t(y)=S_t(y)\sqrt{p_t(y)}/2$. Along a sufficiently small fixed-direction path, $\partial_t F_t/F'_t$ is bounded and its $z$ derivative is polynomially bounded. For each fixed $y>0$, $f_t(y)\to f_0(y)$. Moreover $\|f_t\|_2^2=\frac14\int S_t(F_t(z))^2\phi(z)dz$ converges to $\|f_0\|_2^2$ by a uniform polynomial Gaussian envelope. Pointwise convergence and bounded $L^2$ norms imply weak convergence (first test bounded functions on compact sets, then use density); convergence of norms upgrades this to strong $L^2$ convergence. Integrating this $L^2$-continuous square-root-density derivative proves DQM. This is a statement for each fixed admissible path, not a uniform nonlinear expansion on an $X$-ball.

## 2. Target derivative, independently separated

Let $f_h(x)=\phi(h^{-1}x)/h'(h^{-1}x)$, and define

\[
K_1(x)=\phi(x-\alpha_1-\alpha_2)-\phi(x-\alpha_1)-\phi(x-\alpha_2)+\phi(x),
\]
\[
K_2(x)=f_h(x-\beta_1-\beta_2)-f_h(x-\beta_1)-f_h(x-\beta_2)+f_h(x).
\]

Changing variables separately in the four expectations gives

\[
DI[v]=T_g(r)+b_c^Td,\qquad T_g(r)=\int K_1r_1+\int K_2r_2.
\tag{3}
\]

All these integrals are absolutely convergent. Set $m_0=e^{1/2}$, $m_h=E e^{h(Z)}<\infty$. Differentiating the original contrast, rather than holding its kernels fixed while varying coefficients, gives

\[
b_c=\left(
m_0e^{\alpha_1}(e^{\alpha_2}-1),\;
m_0e^{\alpha_2}(e^{\alpha_1}-1),\;
m_he^{\beta_1}(e^{\beta_2}-1),\;
m_he^{\beta_2}(e^{\beta_1}-1)
\right)^T.
\tag{4}
\]

These four terms are not oracle observations; they describe the derivative at the chosen truth.

## 3. A contracting identity resolves the profile obstruction

**PROVED HERE — concrete model-specific step.** Define

\[
T_j(x)=h^{-1}(h(x)+j),\qquad C(x)=T_2(x)-1,\qquad J=[3,6].
\]

For the anchor response perturbations put

\[
E(x)=q_0(T_2x)+q_0(T_1x)+q_0(Cx)
      -q_4(x)-q_4(T_1x)-q_3(Cx).
\tag{5}
\]

Direct substitution cancels every $r_2$ term and all extraneous $r_1$ terms:

\[
E(x)=r_1(Cx)-r_1(x).
\tag{6}
\]

Anchors have no coefficient derivatives, so (6) also holds for mixed directions. This is an identity of responses inferred from their laws, with no observed cross-stratum coupling.

The interval is invariant: $h(4)-h(3)<2<h(7)-h(6)$ implies $3<C(3)<C(6)<6$. Moreover

\[
0<C'(x)\le\lambda:=575/613<1\quad(x\in J).
\tag{7}
\]

Here is an analytic bound, independent of the diagnostic grid. For $t=T_2x\in(3,7)$, $h''>19/100$ on $[3,7]$, $h'<5/2$ there, and $h'(x)<23/10$ for $x\le6$. These follow directly from

\[
h'=1+\tfrac1{10}(1+2x^2)/\sqrt{1+x^2},\quad
h''=\tfrac1{10}x(3+2x^2)/(1+x^2)^{3/2},\quad
h'''=\tfrac3{10}(1+x^2)^{-5/2}.
\]

Thus $t-x\ge4/5$, $h'(t)-h'(x)\ge19/125$, and $C'=h'(x)/h'(t)\le(23/10)/(23/10+19/125)=575/613$.

For any $x,y\in J$, telescoping (6) yields

\[
r_1(x)-r_1(y)=-\sum_{k=0}^{\infty}
   [E(C^kx)-E(C^ky)].
\tag{8}
\]

No fixed-point profile value is presumed known. The terminal difference tends to zero by continuity and contraction. The stronger operator convergence needed for score continuity is proved below.

The derivative equivalent is useful only for deterministic validation. On $L^2(J)$, set $Bf=C'f\circ C$. Changing variables shows $\|B\|\le\sqrt\lambda=:\rho<1$. Differentiating (6) therefore gives

\[
r'_1=-\sum_{k\ge0}B^kE',\qquad
\left\|r'_1+\sum_{k=0}^{K-1}B^kE'\right\|_2
\le\frac{\rho^K}{1-\rho}\|E'\|_2.
\tag{9}
\]

This argument does not require differentiating an empirical quantile.

## 4. Outward telescoping and finite individual source norms

Eliminate $r_2$ using $r_2(hz)=q_0(z)-r_1(z)$. Define

\[
A_*(z)=h'(z)K_2(h(z)),\quad L(z)=K_1(z)-A_*(z),\quad\int L=0,
\]
\[
P_L(t)=\sum_{n\in\mathbb Z}L(t+n),\quad 4\le t\le5,
\qquad \int_4^5P_L=0,
\]
\[
W(z)=\begin{cases}
\sum_{k\ge1}L(z+k),&z\ge4,\\
-\sum_{k\ge0}L(z-k),&z<4.
\end{cases}
\]

Using $q_3-q_0=r_1(z+1)-r_1(z)$, sum integer shifts toward the central cell $[4,5]$. Absolute convergence justifies the interchange and gives the exact identity

\[
T_g=\int[(A_*-W)q_0+Wq_3]\,dz
 -\sum_{k\ge0}\int_4^5 P_L(t)
      [E(C^kt)-E(C^k4)]\,dt.
\tag{10}
\]

To check its sign, for $z=t+n,n>0$, $r_1(t+n)=r_1(t)+\sum_{j=0}^{n-1}(q_3-q_0)(t+j)$, which gives the positive-tail weight $\sum_{k\ge1}L(z+k)$. Negative $n$ gives the minus sign in $W$. The central contribution is $\int_4^5P_L(t)[r_1(t)-r_1(4)]dt$, and (8) gives (10).

Unlike a one-sided infinite anchor formula, $W$ decays at both infinities. The periodic residue is retained on a compact interval and resolved by (8). This is not cancellation between independent source variances.

For a centered source score the quantile derivative is

\[
q_a(z)=-\frac{F'_a(z)}{\phi(z)}\int_{-\infty}^z S_a(t)\phi(t)dt.
\tag{11}
\]

Indeed $(u_a\phi)'=-S_a\phi$, with zero boundary. Define the vector $\kappa_{a,z}\in\mathcal H$, with only component $a$ nonzero, by

\[
\kappa_{a,z}(t)=\frac{F'_a(z)}{\pi_a\phi(z)}
              [\Phi(z)-\mathbf1\{t\le z\}].
\tag{12}
\]

It is centered, represents $q_a(z)$, and

\[
\|\kappa_{a,z}\|_{\mathcal H}^2
=\frac{F'_a(z)^2\Phi(z)(1-\Phi(z))}{\pi_a\phi(z)^2}.
\tag{13}
\]

On each compact interval, $\|\kappa_{a,x}-\kappa_{a,y}\|\le C_K|x-y|^{1/2}$. This follows from bounded smoothness of $F'_a/\phi$ and the variance of a centered indicator increment. Let $e_x$ be the same six-term combination of $\kappa$'s as (5). The maps $T_j,C$ are Lipschitz on $J$, so

\[
\|e_x-e_y\|\le C_E|x-y|^{1/2}.
\tag{14}
\]

An explicit profile influence is the following convergent Bochner integral and series in $\mathcal H$:

\[
\boxed{\displaystyle
\psi_A=\int[(A_*-W)\kappa_{0,z}+W\kappa_{3,z}]\,dz
-\sum_{k\ge0}\int_4^5P_L(t)[e_{C^kt}-e_{C^k4}]\,dt.}
\tag{15}
\]

The core has norm at most

\[
D_{\rm core}=\frac{C_E}{1-\rho}
       \int_4^5|P_L(t)|\sqrt{|t-4|}\,dt<\infty;
\tag{16}
\]

truncating after $K$ terms leaves at most $\rho^K D_{\rm core}$. The tails are also finite individually. Specifically, for some finite constants and integer $m$,

\[
|A_*|+|L|+|W|\le C(1+|z|)^m e^{-z^2/2+C|z|},\quad
F'_a(z)\le C(1+|z|)^m e^{\delta_0 z^2+C|z|}.
\tag{17}
\]

For the first bound, the $K_1$ terms are Gaussian shifts. For a term $h'(z)f_h(h(z)-b)$, write $w=h^{-1}(h(z)-b)$; since $h'\ge1$, $|w-z|\le|b|$, and $h'(z)/h'(w)$ has a polynomial bound. Thus the same Gaussian bound holds. Summing outwards preserves this tail form. Mills bounds in (13) imply

\[
\|\kappa_{a,z}\|\le C(1+|z|)^m F'_a(z)e^{z^2/4}.
\]

Consequently

\[
D_{\rm tail}:=\int\bigl(|A_*-W|\|\kappa_{0,z}\|+
                                    |W|\|\kappa_{3,z}\|\bigr)dz<\infty,
\tag{18}
\]

because the final quadratic exponent is $(\delta_0-1/4)z^2=-.15z^2$. This also supplies explicit tail error bounds by integrating the displayed envelope over omitted regions. Constants can be bounded from the explicit warp and finite shifts; no numerical variance estimate is used to assert finiteness.

It follows that $\|\psi_A\|\le D_{\rm tail}+D_{\rm core}$ and $\langle\psi_A,S_g(r)\rangle=T_g(r)$ for **every** original analytic direction. This equality extends to the specified completion by continuity; the nuisance closure below is the full closure in $\mathcal H$. The construction uses centered step functions, so it does not mistake smooth weight equations for necessary conditions, nor exclude normalization-related annihilators.

All five coordinates have been specified: $\psi_A$'s coordinates 1 and 2 are zero. This is a stronger anchor sufficiency result for the profile subproblem, not a replacement observation model. The proposed extra laws need not cancel the periodic tail: (5) already removes it using the nonlinear warp. Laws 1 and 2 will supply the four coefficient constraints. Every coordinate of the final influence has finite $L^2$ norm; independence has nowhere been used to cancel a variance.

## 5. Constructive coefficient duals and full nuisance closure

First construct representers $R_1(x)$ of $r_1(x)$ using anchors only. Normalization gives

\[
r_1(4)=\sum_{j=0}^3[q_3(j)-q_0(j)].
\tag{19}
\]

For $x\in[4,5]$, combine (19) and (8), replacing $q$ by $\kappa$ and $E$ by $e$. The resulting $R_1(x)\in\mathcal H$ is a finite sum plus a convergent geometric series. For any other finite $x$, choose an integer $n$ with $x+n\in[4,5]$ and transport back using the finite identity $r_1(z+1)-r_1(z)=q_3(z)-q_0(z)$. Then set

\[
R_2(y)=\kappa_{0,h^{-1}y}-R_1(h^{-1}y).
\tag{20}
\]

These are bounded score functionals representing both profile evaluations even for mixed directions. In particular they pair to zero with each pure coefficient score $s_i$, since they have only anchor components.

For $i=1,2$, let $h_1=h(1)=1+\sqrt2/10$, and form two vector-valued residuals

\[
V_{i,\pm}=\kappa_{i,\pm1}-R_1(\pm1+\alpha_i)-R_2(\pm h_1+\beta_i).
\]

Their pairing with a full score is

\[
\begin{pmatrix}\langle V_{i,-},S(v)\rangle\\
\langle V_{i,+},S(v)\rangle\end{pmatrix}
=M_i\begin{pmatrix}d_{\alpha_i}\\d_{\beta_i}\end{pmatrix},\quad
M_i=\begin{pmatrix}
e^{-1+\alpha_i}&e^{-h_1+\beta_i}\\
e^{1+\alpha_i}&e^{h_1+\beta_i}
\end{pmatrix}.
\tag{21}
\]

The exact determinant is

\[
\det M_i=2e^{\alpha_i+\beta_i}\sinh(\sqrt2/10)>0.
\tag{22}
\]

Applying $M_i^{-1}$ to the two $V$'s defines $\chi_{\alpha_i},\chi_{\beta_i}\in\mathcal H$ satisfying

\[
\langle\chi_j,S_g(r)\rangle=0\ \text{for all }r,
\qquad \langle\chi_j,s_i\rangle=\delta_{ji}.
\tag{23}
\]

**PROVED HERE — the compatibility issue is resolved on the full closure.** Let $\mathcal H_0=(\overline{\operatorname{Range}S_g})^\perp$, $q_i^{\rm eff}=P_{\mathcal H_0}s_i$, and $J_{\rm eff}[i,j]=\langle q_i^{\rm eff},q_j^{\rm eff}\rangle$. Equation (23) and continuity put every $\chi_j$ in $\mathcal H_0$. Let $\Lambda f=(\langle\chi_j,f\rangle)_j$. For every $d$,

\[
|d|=\left|\Lambda\sum_i d_i q_i^{\rm eff}\right|
\le\|\Lambda\|\left\|\sum_i d_i q_i^{\rm eff}\right\|,
\qquad J_{\rm eff}\succeq\|\Lambda\|^{-2}I_4>0.
\tag{24}
\]

No finite profile-bank projection or numerical rank assertion is needed. The bound is finite and constructive, not asserted to be well conditioned in practice.

Because $\psi_A$ has only anchor coordinates, $\langle\psi_A,s_i\rangle=0$. A fully explicit, not necessarily minimum-norm solution is

\[
\boxed{\psi=\psi_A+\sum_{j=1}^4 b_{c,j}\chi_j\in\mathcal H,
\quad DI[v]=\langle\psi,S(v)\rangle\quad\text{for every admissible }v.}
\tag{25}
\]

Its norm is at most $D_{\rm tail}+D_{\rm core}+\sum_j|b_{c,j}|\|\chi_j\|$. This is an equivalent constructive decomposition of the required score equation, not a formal pseudoinverse restatement.

For comparison, take the minimum profile solution $\psi_g=P_{\overline{\operatorname{Range}S_g}}\psi_A$, $d_{\rm res,i}=b_{c,i}-\langle\psi_g,s_i\rangle$. Equation (24) proves compatibility for every residual vector. **STANDARD RESULT APPLIED:** orthogonal decomposition yields

\[
\psi_{\rm eff}=\psi_g+\sum_iq_i^{\rm eff}(J_{\rm eff}^{-1}d_{\rm res})_i,
\quad\|\psi_{\rm eff}\|^2=\|\psi_g\|^2+d_{\rm res}^TJ_{\rm eff}^{-1}d_{\rm res}.
\tag{26}
\]

The second term is orthogonal to the first. Within $\mathcal H_0$, projection onto the span of the four $q_i^{\rm eff}$ gives the smallest norm satisfying its constraints, proving minimality. We have proved the needed properties of these objects via explicit duals; their numerical entries and the efficient variance have not been computed.

## 6. A controlled, infinite-dimensional neighborhood

**PROVED HERE — bounded local extension, not a global rate.** Keep $h,\pi$ fixed and intersect the original model with

\[
\max_{j\in\{1,2\},\,k\in\{0,1,2,3\}}\sup_x e^{-x}|g_j^{(k)}(x)-e^x|\le\epsilon,
\qquad \max_i|c_i-c_{i,0}|\le\epsilon,\quad0\le\epsilon\le.01.
\tag{27}
\]

This is an explicitly controlled neighborhood, not an open ball in the weaker $X$-norm. Arbitrary sufficiently small original analytic directions enter it. The contraction, anchor reconstruction, and normalization identities do not change. Replace $F'$, target weights, and coefficient target derivatives by their values at the new truth. The envelopes and constants in (14), (17), and (18) are uniform over (27).

At the new truth replace $M_i$'s entries by $g'_1(\pm1+\alpha_i)$, $g'_2(\pm h_1+\beta_i)$. Each entry divided by its reference exponential lies between $e^{-3\epsilon}$ and $e^{3\epsilon}$. Hence

\[
\det M_{i,\eta}\ge
2e^{\alpha_{i,0}+\beta_{i,0}}
\sinh(\sqrt2/10-6\epsilon)>0.
\tag{28}
\]

Thus (25) exists throughout this neighborhood. There is also a useful bound on a particular constructed representer, pulled back to the common Gaussian spaces:

\[
\|\widetilde\psi_\eta-\widetilde\psi_0\|_{\mathcal H}
\le C\sqrt\epsilon,
\tag{29}
\]

for a fixed finite $C$. Details that prevent an unjustified smoothness claim: at fixed evaluation points, the coefficient $F'_\eta/\phi-F'_0/\phi$ and its derivative are $O(\epsilon)$ on every compact interval. Therefore the difference between the two kernels' increments is bounded by $C\epsilon\sqrt{|x-y|}$. Geometrically summing (14) proves an $O(\epsilon)$ operator change for the core reconstruction at fixed arguments. Changing an argument by $O(\epsilon)$ costs $O(\sqrt\epsilon)$; within the central chart this follows directly from $|C^kx-C^ky|\le\lambda^k|x-y|$ in (8). Finite transports and (20) preserve that bound. For the calibration queries, $-1+\alpha_i$ and $h^{-1}(h(-1)+\beta_i)$ stay in $(-1,0)$, while their positive counterparts stay in $(1,2)$. Fixed integer shifts $+5$ and $+3$, respectively, put them in $[4,5]$ throughout (27); no change of chart is needed. The tail weights and kernel coefficients change by $O(\epsilon)$ times the integrable envelope in (18). Finally (28), matrix inversion in dimension two, and the coefficient derivative integrals preserve the bound. This proves (29) for the explicit solution (25), not for the efficient orthogonal projection as a function of truth.

The constants may be very large: some anchor queries are in the Gaussian tail near $z=7$. No useful sample size or small efficient variance follows from finiteness. Also, $\delta_0$ is fixed; the contraction and calibration determinant degenerate toward a different problem as the warp becomes linear. No uniform-in-$\delta$ claim is made.

## 7. Feasibility and the exact outstanding statistical step

At the population level the identities also hold with $g_j,F_a$ in place of $r_j,q_a$. Replace (19)'s base value 0 by $g_1(0)=1$. Thus a concrete possible estimator would use the five empirical quantile functions at $\Phi(z)$, finite versions of the convergent reconstruction, a locally selected solution of the two calibration equations at $z=\pm1$ for each coefficient pair, and then the outward target formula. This description uses observed responses only, not latent draws or target labels. It does not require recovering entire profiles on the line. It has not been implemented or run here.

**OPEN — attainment.** The missing step is a uniform asymptotic expansion for a feasible truncated empirical-quantile/calibration estimator with estimated coefficients:

\[
\widehat I-I(\eta)=\sum_a\pi_a\frac1{n_a}
      \sum_{i=1}^{n_a}\psi_{\eta,a}(Y_{ai})+o_{P_\eta}(N^{-1/2}).
\tag{30}
\]

Finite variance of (25) does not prove (30). A proof must select growing central-series depth $K_N$, tail cutoff, and quadrature accuracy with total deterministic bias $o(N^{-1/2})$; establish stochastic weighted-quantile control for those choices; and control calibration/root selection and the nonlinear remainder uniformly in a suitable neighborhood. Geometric series tails suggest $K_N$ proportional to $\log N$, but that alone controls neither tail empirical quantiles nor the accumulation of their nonlinear errors. Empirical quantiles do not satisfy the analytic profile model. No unproved quantile rate, DKW-to-Hellinger bridge, or generic nuisance convergence rate is assumed.

An alternative one-step method would additionally need feasible estimates of the profiles/derivatives defining $\psi_\eta$, a consistent source-only calibration, an $L^2$-consistent estimated influence under the appropriate pullback, and a bias/remainder bound of order $o_P(N^{-1/2})$; sample splitting could control the empirical term but cannot supply the missing bias bound. Bound (29) describes controlled deterministic perturbations, not a proved nuisance-estimator rate. The influence at the reference remains an oracle analysis object until such estimation is justified.

With the convention $H^2(P,Q)=\int(\sqrt p-\sqrt q)^2$, fixed-path DQM gives $d_\pi(P_{\eta_t},P_\eta)=|t|\|S(v)\|/2+o(|t|)$. Equation (25) therefore bounds each directional target derivative in the complete five-law score metric. It does not supply an endpoint Hellinger Lipschitz bound uniform over pairs, global minimax orders, or a sharp attainable rate. The previous global unresolved statement survives.

## 8. Closest theory and validation scope

**STANDARD RESULT APPLIED.** [Trabs, Theorem 2.7 and its following discussion](https://arxiv.org/html/1307.6610v2) supply the adjoint-range criterion and the efficient influence/information-bound framework for regular indirect experiments; solving that criterion does not itself construct an attaining estimator in the nonlinear experiment. The orthogonal nuisance projection and finite coefficient correction in (26) are the same established machinery. Quantile derivatives, telescoping, and a contraction-series inverse are also standard tools.

The concrete argument supplied here is (5)–(18): the nonlinear known warp creates a contracting private-anchor identity, so outward target weights plus a compact geometric remainder yield a bounded profile functional. Equations (19)–(24) then prove full coefficient compatibility by explicit bounded profile-evaluation and coefficient duals. This is the model-specific solution that an abstract adjoint criterion alone did not give. A claim that this argument is new throughout the literature would need a separate comparison beyond the cached primary sources; no such claim is made. It is a local statistical theorem, not a completed learned-model or ICLR contribution, nor a biological siRNA result.

**NUMERICALLY CHECKED, not existence certificates.** [checks.py](checks.py) and [checks.json](checks.json) validate the six-term identity, derivative contraction sums, point reconstruction/coefficient residuals, and the direct profile functional on two original-class analytic directions, four coefficient basis directions, and one mixed direction. No finite dictionary is used as the tangent model. There is no score Gram solve, spectral variance estimate, sample, fit, or pilot rerun. The analytic proof above supplies the infinite-dimensional result.

The first attempt failed an absolute $2\times10^{-10}$ assertion for a six-term cancellation with terms of order $2\times10^5$. Its script and failure record are preserved. The second attempt prints absolute and scaled residuals and checks the scaled error, without changing a direction, reference, or mathematical formula. At 512 contraction terms, maximum derivative errors are $2.881\times10^{-9}$ and $8.801\times10^{-10}$; the maximum coefficient recovery error is $1.754\times10^{-9}$. At quadrature order 128, the two target identity residuals are $6.484\times10^{-14}$ and $4.797\times10^{-14}$. The 64-to-128 source-integral change is $6.220\times10^{-8}$. These are floating-point diagnostics, not certified quadrature enclosures or upper bounds on the efficient variance. Truncation and both tails are bounded analytically above, separately from solver/roundoff residuals. See [commands and resources](commands_and_resources.md) for the failed attempt, successful job, and allowance charges.
