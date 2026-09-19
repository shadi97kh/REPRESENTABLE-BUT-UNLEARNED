# Certified piecewise integration and arithmetic allocation

This is a separate static derivation for the compact construction. It does not modify the frozen local-attainment result. No numerical job, sampled input, fit, or benchmark was used in deriving this note.

**PROVED HERE, with explicit arithmetic inputs below:** the finite compact-anchor target can be integrated to total absolute error \(O(N^{-1})\), with \(N\,\operatorname{polylog}N\) arithmetic operations. A deterministic implementation of the finite coefficient iteration can be included using \(O((\log N)^2)\) working bits under the stated rounded-observation and certified-primitive model. This is an algorithmic existence/work bound, not certification of the accompanying mpmath diagnostic implementation or evidence of useful constants.

Every occurrence of \(N^{-1}/8\) here means \(1/(8N)\).

## 1. Finite problem and polygon structure

Let \(N=5n\ge10\), \(0<Y_{ai}\le N^2\), and let \(Q_a\) be the polygon through \((i/n,Y_{a(i)})\), constant below \(1/n\). Then
\[
0<Q_a\le N^2,\qquad
\operatorname{Lip}(Q_a)\le nN^2\le N^3.
\]
Inside each polygon segment, \(Q_a(p)=a_a+b_ap\) with \(0\le b_a\le N^3\). All source groups have the same *probability* knots; their observed heights are different.

Use \(J=[-7/4,-3/4]\), \(x_*=-1\), \(\lambda=271/290\),
\[
C(x)=h^{-1}(h(x)+5)-4,\quad
K=\left\lceil2\log N/|\log\lambda|\right\rceil+2,\quad
R=\sqrt{1.5\log N},\quad M=\lceil4R+30\rceil+2 .
\]
The finite target consists of integrals of
\[
U_{0,c,M}(z)Q_0(\Phi z)+W_{c,M}(z)Q_3(\Phi z)
\]
over \([-R,R]\), and the compact integral
\[
P_{c,M}(t)D_K(Q;t,x_*),\qquad t\in J.
\]
Any explicit finite-periodization normalization correction is retained; its integrand has the same type of bounds, including its nonzero finite integral. The known outward-weight jump is split at the unit-cell boundary used in the construction, not at the normalization point merely because that point is \(-1\).

The simplified \(E_{45}\) has 18 signed response terms:
\[
\sum_{j=1}^{5}F_0(T_jx)
+\sum_{r=1}^{4}F_0(T_5x-r)
-\sum_{j=0}^{4}F_4(T_jx)
-\sum_{r=1}^{4}F_3(T_5x-r).
\]
Thus \(D_K=-\sum_{k<K}\{E_{45}(C^kt)-E_{45}(C^kx_*)\}\) contains at most \(36K\) response values. The reference half is constant in the integration variable and can be cached. This is algebra on the five observed quantile curves; transformed arguments are not extra actions.

## 2. Fully public fourth-derivative envelopes

**STANDARD RESULT APPLIED:** write \({\cal B}_{r,m}\) for the partial exponential Bell polynomial in the ordinary derivative convention:
\[
(f\circ g)^{(r)}=\sum_{m=1}^{r}f^{(m)}(g){\cal B}_{r,m}
(g',\ldots,g^{(r-m+1)}).
\]
All Bell coefficients are nonnegative.

Useful global warp bounds are
\[
H_1(B)=1.1+.2B,\quad H_2=.2,\quad H_3=.3,\quad
H_4=1.5,\quad H_5=10.5.
\]
They follow directly from
\[
\begin{split}
h''(x)&=.1x(3+2x^2)(1+x^2)^{-3/2},\\
h'''(x)&=.3(1+x^2)^{-5/2},\\
h''''(x)&=-1.5x(1+x^2)^{-7/2},\\
h'''''(x)&=(-1.5+9x^2)(1+x^2)^{-9/2}.
\end{split}
\]
The inverse \(w=h^{-1}\), using \(h'\ge1\), has global bounds
\[
a_1=1,\quad a_2=.2,\quad a_3=.42,\quad
a_4=2.22,\quad a_5=17.328.
\]
For example these follow inductively from
\[
a_r=\sum_{m=2}^{r}H_m{\cal B}_{r,m}(a_1,\ldots,a_{r-m+1}),
\]
which is a bound on the differentiated inverse identity, not an equality of derivatives.

Put
\[
b_r=\sum_{m=1}^{r}a_m{\cal B}_{r,m}
(H_1(7/4),H_2,\ldots,H_{r-m+1}),\qquad r\le4.
\]
These bound the derivatives of every \(T_j\) on \(J\), and the derivatives of \(C\) of order at least two. Use the sharper proved contraction bound \(|C'|\le\lambda\) at order one. Define
\[
D_1=1,\qquad D_r=(1-\lambda)^{-1}
 \sum_{m=2}^{r}b_m{\cal B}_{r,m}
 (D_1,\ldots,D_{r-m+1}),\quad 2\le r\le4.
\]
Induction in \(k\), followed by induction in \(r\), proves
\[
\sup_{k\ge0,t\in J}|(C^k)^{(r)}(t)|\le D_r.
\]
Consequently all outer queries \(f(t)=T_j(C^kt)-r\), with the shifts actually present in \(E_{45}\), have derivative bounds
\[
f_r=\sum_{m=1}^{r}b_m{\cal B}_{r,m}
(D_1,\ldots,D_{r-m+1}).
\]
The identity outer map is covered by these upper bounds.

Gaussian derivative bounds can be taken as \(d_0=d_1=1\) and
\(d_r=(r-1)d_{r-2}\). The Fourier integral of the Gaussian and the recurrence for absolute Gaussian moments prove \(|\phi^{(r)}|\le d_r\). It follows that
\[
p_r=\sum_{m=1}^{r}d_{m-1}{\cal B}_{r,m}(f_1,\ldots,f_{r-m+1})
\]
bounds the derivatives of \(\Phi(f(t))\). On a polygon piece,
\[
|D_K|\le36KN^2,\qquad
|D_K^{(r)}|\le18KN^3p_r,\quad1\le r\le4.
\]
No derivative of an unknown profile or empirical density enters these bounds.

Here are equally explicit bounds for the weights. Let \(\rho(y)=\phi(w(y))w'(y)\). Set
\[
G_0=d_0,\qquad
G_\ell=\sum_{m=1}^{\ell}d_m{\cal B}_{\ell,m}(a_1,\ldots,a_{\ell-m+1}),
\]
and
\[
{\cal R}_r=\sum_{\ell=0}^{r}{r\choose\ell}G_\ell a_{r-\ell+1},
\quad0\le r\le4.
\]
Then \(|\rho^{(r)}|\le{\cal R}_r\). With
\[
C_0(B)={\cal R}_0,\qquad
C_\ell(B)=\sum_{m=1}^{\ell}{\cal R}_m
 {\cal B}_{\ell,m}(H_1(B),H_2,\ldots),
\]
define
\[
A_r(B)=4\sum_{\ell=0}^{r}{r\choose\ell}
 H_{r-\ell+1}(B)C_\ell(B),\qquad
L_r(B)=4d_r+A_r(B).
\]
For \(r\le4\), these bound \(A^{(r)}\) and \(L^{(r)}\), uniformly over the public coefficient rectangles; translated density derivatives have the same global bounds. Take \(B=R+M+4\). Then
\[
P_r=(2M+1)L_r(B),\quad W_r=ML_r(B),\quad
U_r=A_r(R)+ML_r(B)
\]
bound the finite weights on each side of their jump.

A concrete central fourth-derivative bound is
\[
{\cal C}_{4,\mathrm{cen}}=
36KP_4+18K\sum_{\ell=1}^{4}{4\choose\ell}P_{4-\ell}p_\ell.
\]
A concrete outer bound is
\[
{\cal C}_{4,\mathrm{out}}=
U_4+W_4+\sum_{\ell=1}^{4}{4\choose\ell}
 (U_{4-\ell}+W_{4-\ell})d_{\ell-1}.
\]
If an additional normalization integrand is present, add its corresponding Leibniz bound rather than dropping it. Therefore \(|f^{(4)}|\le N^3{\cal C}_4(N)\) on each smooth piece, with a completely public \({\cal C}_4\) that is a polynomial in \(K,R,M\), with possibly large fixed constants.

Similarly let \(B_0=N^2{\cal C}_0(N)\) bound the absolute integrands, using \(36KP_0\) centrally and \(U_0+W_0\) outside, and including the normalization term if present.

## 3. Locating knots without assuming stable recursive inverses

Every source query in the central terms lies in \([-7/4,13/4]\). On this range \(1<h'<2\), hence
\[
(T_j)'>1/2,\qquad (C^k)'\ge2^{-k}.
\]
Thus every varying central probability map \(p_{jk}(t)\) has
\[
p_{jk}'(t)\ge m_K:=2^{-K-11}.
\]
Indeed \(\phi(13/4)>2^{-10}\), so this lower bound is conservative. Conversely,
\(p_{jk}'(t)\le2\lambda^k\). A monotone probability map of range length at most \(2\lambda^k\) crosses at most \(2n\lambda^k+2\) knots. Summing over 18 terms and all levels, and adding the outer probability knots, gives the public overcount
\[
B_{\mathrm{knots}}\le
36n/(1-\lambda)+36K+2n+O(1).
\]
In particular the total is \(O(N+K)\), not \(O(NK)\). Duplicate preimages may be merged but need not be identified as exactly equal.

Enumerate candidate ranks from certified endpoint probability enclosures. Endpoint uncertainty can include neighboring ranks; it does not exclude a possible root. To localize a root \(p_{jk}(t)=i/n\), use forward bisection in \(t\), not repeated evaluation of \(C^{-1}\). For desired enclosure width \(\delta\), evaluate the probability to error at most \(m_K\delta/16\). If the result is separated from \(i/n\), choose the corresponding half interval. If it straddles \(i/n\), the derivative lower bound places the root inside an interval of width at most \(\delta/4\) around the midpoint. Thus exact equality cannot cause nontermination. The number of bisection steps is \(O(\log(1/\delta))\).

For outer knots, the lower derivative bound is \(\phi(R)\) on \([-R,R]\); it is \(N^{-3/4}/\sqrt{2\pi}\), also polynomially small. Known jump points are inserted directly. Endpoint and node enclosures are treated using the same error allocation.

## 4. Gaps and certified Gauss quadrature

Allocate \(\varepsilon_{\mathrm{gap}}=\varepsilon_{\mathrm{quad}}=1/(8N)\).
Choose a dyadic enclosure width
\[
\delta\le
\frac{\varepsilon_{\mathrm{gap}}}
 {4(1+B_{\mathrm{knots}})B_0}.
\]
Take the union of all knot enclosures and omit these small gaps from integration. The absolute error is bounded by \(B_0\) times their total length, hence at most \(\varepsilon_{\mathrm{gap}}\), including allocated endpoint enclosures. This step prevents an unproved smoothness assertion at an approximately located knot. It also handles arbitrarily close or overlapping knots.

On the complementary intervals all polygon branches are fixed and the integrand is \(C^4\). Split further into panels of length at most
\[
h_N=\min\left\{1,
\left(\frac{4320\varepsilon_{\mathrm{quad}}}
 {N^3{\cal C}_4(N)(1+2R)}\right)^{1/4}\right\}.
\]
The two-node Gauss rule has error at most
\(\|f^{(4)}\|_\infty\ell^5/4320\) on an interval of length \(\ell\). Summation proves total quadrature error at most \(\varepsilon_{\mathrm{quad}}\).

The panel count is
\[
O\!\left(B_{\mathrm{knots}}+(1+2R)/h_N\right)
=N\,\operatorname{polylog}N .
\]
A central evaluation costs \(O(K+M)\) quantile/inverse operations when the iterates are propagated once and the fixed reference contribution is cached. Knot localization costs \(O(K\log N)\) primitives per central root using the direct forward maps. Sorting, preimage localization, panel integration and weight evaluation consequently use \(N\,\operatorname{polylog}N\) primitives. Memory is \(O(N+K+M)\) words with panel streaming; retaining all panels unnecessarily increases memory.

**Important limitation:** the displayed constants were chosen for easy certified bounds. A near-linear asymptotic work count is not evidence that these conservative constants are affordable.

## 5. Coefficient iteration and a complete finite arithmetic comparison

This section compares two *finite* statistics: the compact target using the exact finite projected iteration, and its rounded arithmetic implementation. Whether that exact iteration has the desired statistical expansion is proved separately; a numerical comparison cannot supply a missing statistical iteration argument.

Suppose the exact iteration uses the public preconditioner for \(T=O(\log N)\) steps. On every fixed calibration chart verified over the complete public rectangles, the capped polygon reconstruction has a public bound
\[
\operatorname{Lip}(\widehat B)\le C_BN^3(K+1).
\]
Finite translations and the known inverse warp only change \(C_B\). The constant is obtained from the finitely many query derivatives and the absolute counts of their signed coefficients. Projection is nonexpansive. Therefore
\[
L_N=1+\|M_0^{-1}\|C_BN^3(K+1),\qquad
A_T=\sum_{s=0}^{T-1}\max(1,L_N)^s
\]
bounds accumulation of arbitrary deterministic per-step perturbations. This does **not** assume empirical contraction.

Let \(J_c\) bound the coefficient derivative of the finite target on the public rectangles:
\[
J_c\le N^2 C_c(N),
\]
where \(C_c\) is a public polynomial in \(K,M,R\), obtained by differentiating the known translated kernels. Put
\[
\delta_c=\frac{1}{16N(1+J_c)}\quad\text{for each of the two pairs}.
\]
If the total certified error of each projected update is at most
\(\delta_c/A_T\), each final coefficient pair differs from its exact finite-iteration counterpart by at most \(\delta_c\); the complete four-vector error is at most \(2\delta_c\); the resulting target error is at most \(1/(8N)\). The update budget includes errors in the public matrix inverse, each source response, normalization, finite reconstruction, summation, and projection inputs.

Here is an explicit rounded-observation allocation. Let \(C_{\rm data}\) be the public sup-norm sensitivity of one update to perturbing every observed height by at most one. The polygon order-statistic operator is 1-Lipschitz in the sup norm, even if rounding changes the ordering. Finite reconstruction gives \(C_{\rm data}\le C(K+1)\). Round each observed response to error
\[
\tau_Y\le
\min\left\{\frac{\delta_c}{4A_TC_{\rm data}},
\frac{1}{32N(1+J_Y)}\right\},
\]
where \(J_Y\le C(K+1)(M+1)(R+1)^m\) is a public height sensitivity of the finite target at fixed coefficients. Reserve the other three quarters of each update error for preconditioner, primitives, and arithmetic. This is an observed-input precision requirement, not access to a latent quantity.

A fixed exponent \(N^{-D}\) is not justified by this worst-case iteration calculation: \(A_T\) can be \(\exp[O((\log N)^2)]\). Instead
\[
\log(1/\tau_Y)=O((\log N)^2).
\]
The theorem must state this input model. On the response-cap event, use upward dyadic rounding followed by clipping at \(N^2\); this keeps every encoded observation strictly positive and at most the cap, with the same absolute error bound. A coarser externally fixed observation precision does not satisfy the comparison proved here. If cap classification itself uses approximate rather than exact input, allocate a guard band and account separately for its probability; do not silently assume an exact real comparison.

All other arithmetic allocations are finite and explicit:

- use probability enclosures of radius at most \(m_K\delta/16\) for central knot decisions, and \(\phi(R)\delta/16\) for outer decisions;
- give integration node and weight errors combined absolute budget \(1/(8N)\);
- give final summation and endpoint errors combined absolute budget \(1/(8N)\);
- bound each straight-line evaluation's error by a certified interval enclosure whose width is below its assigned share, using the public slope \(N^3\), the envelopes above, and the known numbers of terms;
- allocate input and coefficient shares as above. Unused shares remain slack; they are not used to excuse any unchecked operation.

The finite computational graph has no unknown conditioning constant. Bounds for every node can also be computed recursively: for sums add sensitivity bounds; for a product use \(|u|\,d_v+|v|\,d_u\); for a known smooth primitive multiply by its public derivative bound. Use continuous polygon interval evaluation at uncertain ranks, rather than an exact transcendental comparison. Inverse \(h\) evaluation is conditioned by \(h'\ge1\), and contraction iteration query errors sum geometrically. The sole larger amplification is the explicitly retained \(A_T\).

## 6. Bit complexity and scope

**PROVED HERE in the stated arithmetic model:** all arguments to known primitives have magnitude polynomial in \(\log N\), and all intermediate magnitudes have logarithms polynomial in \(\log N\). The required absolute precision is \(O((\log N)^2)\) bits for the complete algorithm above; target integration alone uses \(O(\log N)\). Rational bisection supplies square roots and the known monotone inverse. Range-reduced Taylor series with explicit tails supply the exponential and Gaussian integral in work polynomial in the requested bit count and argument range. Rational rounding after operations prevents denominator lengths growing with the number of summands. A rational series for \(\pi\), or any separately certified elementary-constant routine, completes the primitive model. Thus ordinary polynomial-time integer arithmetic gives a total \(N\,\operatorname{polylog}N\) bit-operation bound. Reading higher-precision input is included; no unit-cost arbitrary-real observation oracle is being claimed.

This is a specified certifiable algorithm, not a statement that mpmath returns enclosures. The constants and enclosure routines must be instantiated before claiming an implemented certified runtime or numerical certificate.

Adding the allocated errors gives a deterministic \(O(N^{-1})\) difference from the same compact finite-iteration ideal statistic on the response-cap event, hence
\[
\sqrt N\,|\widehat I_{\mathrm{arith}}-
\widehat I_{\mathrm{compact,\ finite\ iteration}}|\longrightarrow0 .
\]
The coefficient iteration's statistical equivalence, central truncation, finite periodization remainder, response tails, and local-alternative expansion are separate obligations. This comparison does not assert equality with the old anchor's influence.

The central probability range has minimum tail probability
\(\Phi(-13/4)\), rather than the old \(\Phi(-7)\). The response-tail condition still involves \(\Phi(-(R+1))\). Both are sufficient proof conditions, not operational zero-output gates or necessary sample-size bounds. No finite-sample accuracy or practical validated regime follows from the operation bound.

