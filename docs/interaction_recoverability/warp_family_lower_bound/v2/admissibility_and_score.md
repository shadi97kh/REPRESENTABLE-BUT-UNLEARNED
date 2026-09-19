# Global admissibility, actual score and target separation

> Version 2: the supplied seed and the existing v1 derivation were checked algebraically, including global transition bounds, the common-coordinate density derivative and the target margin. The proof is retained with explicit verification details added below. This is one-agent static self-review, without numerical evaluation or an external review certificate.

Every bound below holds uniformly for \(\theta\in[1/2,2]\), \(0<\delta\le1/3{,}200{,}000\), and \(|t|\le1/20\), unless a fixed-parameter qualification is stated. Notation is defined in [construction.md](construction.md). This is an analytic derivation; no numerical evaluations or diagnostics were run.

## 1. Saturation and the global derivative bounds

**PROVED HERE.** Direct differentiation gives
\[
L_\delta'(v)=\tfrac12\left[
\tanh\frac{v+A}{2\delta}-\tanh\frac{v-A}{2\delta}\right],
\]
\[
L_\delta''(v)=\frac1{4\delta}\left[
\operatorname{sech}^2\frac{v+A}{2\delta}
-\operatorname{sech}^2\frac{v-A}{2\delta}\right].
\tag{1}
\]
The function is odd and increasing with limits \(\pm A\). Put
\(\chi(v)=\exp[-(|v|-A)_+/\delta]\). Using the hyperbolic-function expressions, and \(1-\tanh x\le2e^{-2x}\), \(\operatorname{sech}^2x\le4e^{-2|x|}\), gives
\[
|L_\delta(v)|\le A,\qquad 0<L_\delta'(v)<1,\qquad
L_\delta'(v)\le\chi(v),\qquad
|L_\delta''(v)|\le\delta^{-1}\chi(v).
\tag{2}
\]
For \(|v|\le A\), an exact form is
\[
L_\delta(v)=v+\delta\left[
\log(1+e^{-(A+v)/\delta})
-\log(1+e^{-(A-v)/\delta})\right].
\]
Therefore, with \(E_\delta=e^{-A/(2\delta)}\), on \(|v|\le A/2\),
\[
|L_\delta(v)-v|\le2\delta E_\delta,\qquad
|1-L_\delta'(v)|\le2E_\delta.
\tag{3}
\]
The candidate constants are valid, though deliberately loose.

The warp satisfies
\[
s_\theta'(z)^2=4\theta z^2+\frac1{1+\theta z^2}
\le1+4\sqrt\theta\,|s_\theta(z)|,\qquad
|s_\theta''(z)|\le2\sqrt\theta.
\tag{4}
\]
For the second inequality, writing \(r=\theta z^2\), its squared comparison is
\(r(3+2r)^2\le4(1+r)^3\); the difference is \(4+3r>0\).
Hence
\[
v'(z)^2\le\delta^2+4\delta\sqrt\theta\,|v(z)|,\qquad
|v''(z)|\le2\delta\sqrt\theta.
\]
Since \(\delta\le A\), the functions
\(r e^{-(r-A)_+/\delta}\) and \(r e^{-2(r-A)_+/\delta}\) have supremum at most \(A\). This follows by differentiating on \(r>A\); they decrease there.
Applying the chain rule \(u'=L'(v)v'\), \(u''=L''(v)v'^2+L'(v)v''\), gives the global bounds
\[
0\le u'(z)\le
\sqrt{\delta^2+4\delta\sqrt\theta A}<1/1000,
\]
\[
|u''(z)|\le\delta+4\sqrt\theta A+2\delta\sqrt\theta
\le3/5+4\delta<2/3.
\tag{5}
\]
For the final rational comparisons use \(\sqrt\theta\le\sqrt2<3/2\) and \(\delta<10^{-6}\).
These are global bounds, including the moving transition region. No complex-analytic derivative estimate is assumed.

The normalizing integral in construction (1) lies in \([e^{-A},e^A]\), so \(e^{-A}\le b_\delta\le e^A\). The integrand is positive and real analytic at every real point: real \(\cosh\) is strictly positive, so its logarithm is real analytic, as are the square root and the compositions involved here. Its primitive is locally analytic and finite from \(-\infty\), by comparison with \(e^y\).
Consequently \(g_{1,\delta}\) is positive and real analytic on the entire real line, \(g_{1,\delta}(0)=1\), and \(g_{1,\delta}(-\infty)=0\). Real analyticity on \(\mathbb R\) does not mean complex-entire analyticity.

Direct derivatives are
\[
g_1'(x)=b_\delta e^{x+u(x-a_0)},\quad
g_1''(x)=g_1'(x)[1+u'(x-a_0)],
\]
\[
g_1'''(x)=g_1'(x)\{[1+u'(x-a_0)]^2+u''(x-a_0)\}.
\tag{6}
\]
Using \(e^{1/5}<5/4\), which follows from the geometric upper bound on its positive series,
\[
\tfrac45e^x<g_1'(x)<\tfrac54e^x,\quad
|g_1''(x)|\le\tfrac54(1001/1000)e^x<2e^x,
\]
\[
|g_1'''(x)|\le\tfrac54\{(1001/1000)^2+2/3\}e^x<3e^x.
\tag{7}
\]
Thus every original derivative envelope holds with strict margins. The second profile \(e^x\) plainly satisfies the same original class.

## 2. Coefficients and the exact response derivative

**PROVED HERE.** For \(|t|\le1/20\), \(b_\delta\le e^{.1}\le10/9\),
\(e^{.05}-1\le1/19\), and \(1-e^{-.05}\le1/20\). These imply
\[
161/171\le D(t)\le19/18.
\]
Moreover
\(-\log(161/171)\le10/161<1/16\) and
\(\log(19/18)\le1/18<1/16\).
This proves the coefficient bounds in construction (5).
Define
\[
\lambda_t=b_\delta e^{a_0+t-\beta_1(t)}
=b_\delta e^t/D(t).
\]
Elementary lower bounds \(e^{-.1}\ge9/10\) and \(e^{-.05}\ge19/20\) also give
\[
81/100\le\lambda_t\le200/161<5/4<2,\qquad
\beta_1'(t)=-\lambda_t.
\tag{8}
\]

With all profiles fixed along the path, differentiation gives exactly
\[
\dot F_{1,t}(z)=b_\delta e^{z+a_0+t}
\{e^{u(z+t)}-e^{v(z)}\}.
\]
Put \(d_t(z)=v(z)-u(z+t)\), \(H(z)=h'(z)\). Dividing by
\[
F_{1,t}'(z)=b_\delta e^{z+a_0+t+u(z+t)}
+H(z)e^{z+v(z)+\beta_1(t)}
\]
proves
\[
w_t(z)=\frac{\dot F_{1,t}(z)}{F_{1,t}'(z)}
=\frac{\lambda_t(1-e^{d_t(z)})}
{\lambda_t+H(z)e^{d_t(z)}}.
\tag{9}
\]
This formula preserves the cancellation from the exact nonlinear amplitude balance.

The curve is a strictly increasing smooth bijection \(\mathbb R\to(0,\infty)\), since its derivative is positive, both summands vanish at negative infinity, and \(g_1(z+\alpha_1)\) alone diverges at positive infinity. Its density is strictly positive on that common support:
\[
p_t(y)=\frac{\phi(z_t(y))}{F_{1,t}'(z_t(y))},
\qquad z_t(y)=F_{1,t}^{-1}(y).
\]
The fixed-outcome score follows by differentiating this density and \(\partial_t z_t=-w_t(z_t)\):
\[
\boxed{S_t(z)=z\,w_t(z)-w_t'(z).}
\tag{10}
\]
The Jacobian term \(-w_t'\) is retained. No response-distance or Wasserstein implication is used.

## 3. Central score bound

**PROVED HERE.** Set
\[
r_\delta^2=A/(16\delta),\qquad
c_0=A/128=1/1280,\qquad
\epsilon_\delta=e^{-c_0/\delta}.
\tag{11}
\]
The elementary bounds
\[
|s_\theta(z)|\le |z|+\tfrac32z^2,\quad
s_\theta'(z)\le1+3|z|,\quad |s_\theta''(z)|\le3
\tag{12}
\]
are uniform in \(\theta\). Since \(r_\delta\ge1\) and \(|t|\le1/20\), on \(|z|\le r_\delta\) one has \(|z+t|\le2r_\delta\). Also \(\delta\le A/16\); therefore
\[
|v(z+t)|\le2\delta r_\delta+6\delta r_\delta^2
\le A/8+3A/8=A/2.
\tag{13}
\]
The same holds at \(z\). Equations (3), (12) and the mean-value theorem now give
\[
|d_t(z)|\le4\delta|t|(1+|z|)+2\delta E_\delta,
\]
\[
|d_t'(z)|\le3\delta|t|+8\delta E_\delta(1+|z|).
\tag{14}
\]
In the derivative bound write
\(v'(z)-L'(v(z+t))v'(z+t)\) as its shift difference plus the residual \(1-L'\). This retains the small central residual.

For fixed \(\lambda=\lambda_t>0\), differentiation of
\(w(d,H)=\lambda(1-e^d)/(\lambda+He^d)\) gives
\[
|\partial_dw|
=\frac{\lambda(\lambda+H)e^d}{(\lambda+He^d)^2}
\le\frac{\lambda+H}{4H}\le3/4,
\]
\[
|\partial_Hw|
=\frac{e^d}{\lambda+He^d}|w|
\le |w|/H.
\tag{15}
\]
The first bound follows by maximizing \(x/(1+x)^2\). In the second, cancellation in \(w\) is essential. Since \(w(0,H)=0\), \(|w|\le(3/4)|d|\), and since \(H\ge1\), \(|H'|\le3\delta\),
\[
|w'|\le(3/4)|d'|+3\delta|w|.
\]
Substituting (14), and using \(\delta<1/9\), yields on the central region
\[
|S_t(z)|\le
4\delta|t|(1+|z|)^2+
8\delta E_\delta(1+|z|).
\]
For example the coefficient of \(\delta|t|\) before simplification is
\(3|z|(1+|z|)+9/4+9\delta(1+|z|)\), bounded by the displayed polynomial; the coefficient of \(\delta E_\delta\) is at most \(8(1+|z|)\).
Minkowski and Gaussian moments give
\(\|(1+|Z|)^2\|_2\le3+\sqrt3<5\),
\(\|1+|Z|\|_2\le2\). Consequently the conservative bound
\[
\|S_t\mathbf1_{\{|Z|\le r_\delta\}}\|_2
\le30\delta(|t|+E_\delta)
\tag{16}
\]
holds uniformly.

## 4. Tail score bound and explicit exponent

**PROVED HERE.** Globally, (9) gives \(|w|\le\max(1,\lambda/H)\le2\). From \(0\le L'\le1\) and (12),
\[
|d_t'(z)|\le|v'(z)|+|u'(z+t)|
\le7\delta(1+|z|).
\]
Equations (15) imply \(|w'|\le12\delta(1+|z|)\) and hence
\[
|S_t(z)|\le3(1+|z|).
\tag{17}
\]
For any \(r\ge0\),
\[
E[(1+|Z|)^2\mathbf1_{\{|Z|>r\}}]
\le e^{-r^2/4}E[(1+|Z|)^2e^{Z^2/4}]
<9e^{-r^2/4}.
\]
The last inequality follows from
\((1+|Z|)^2\le2(1+Z^2)\),
\(Ee^{Z^2/4}=\sqrt2\), and
\(EZ^2e^{Z^2/4}=2\sqrt2\).
Thus
\[
\|S_t\mathbf1_{\{|Z|>r_\delta\}}\|_2
\le9e^{-r_\delta^2/8}=9\epsilon_\delta.
\]
Since \(\delta E_\delta\le\epsilon_\delta\), combining with (16) proves
\[
\boxed{\|S_t\|_{L^2(\phi)}
\le40(\epsilon_\delta+\delta|t|).}
\tag{18}
\]
This verifies the requested positive exponent \(c_0=A/128\) with explicit constant 40 independent of \(\delta,\theta,t\). It is not an optimal exponent.

At \(t=0\) the information is small but not zero. Indeed \(0<L'(v)<1\) and \(L(0)=0\) imply \(L(v)<v\) for \(v>0\). Thus \(w_0(z)<0\) for \(z>0\). If \(S_0=0\) almost everywhere, continuity and \((\phi w_0)'=-\phi S_0\) would make \(\phi w_0\) constant. Its tail limit is zero because \(|w_0|\le2\), a contradiction. Consequently \(0<\|S_0\|_2<\infty\). This construction does not assert exact singular information at any fixed positive \(\delta\).

## 5. DQM and control along the entire likelihood path

**PROVED HERE.** To justify a square-root-density integral rather than just a pointwise score calculation, fix \(\delta,\theta\), and write
\[
Z_t(z)=F_{1,t}^{-1}(F_{1,0}(z)),\qquad J_t(z)=\partial_z Z_t(z)>0.
\]
The flow equations are
\[
\partial_t Z_t=-w_t(Z_t),\qquad
\partial_t\log J_t=-w_t'(Z_t),\qquad Z_0=z,\ J_0=1.
\]
By the global bounds above,
\[
|Z_t(z)-z|\le2t_0=1/10,\qquad
J_t(z)\le\exp\{\delta(1+|z|)\}.
\]
For the second inequality integrate \(12\delta(1+|Z_s|)\) over a path of length at most \(1/20\), and use \(12(1/20)(1+|z|+.1)\le1+|z|\).
The Gaussian ratio then gives the convenient uniform bound
\[
\phi(Z_t(z))J_t(z)\le e^{1+|z|}\phi(z).
\tag{19}
\]
The square-root density on a common fixed measure is
\[
f_t(z)=\sqrt{p_t(F_{1,0}(z))F_{1,0}'(z)}
=\sqrt{\phi(Z_t(z))J_t(z)}.
\]
Its derivative is \(\partial_t f_t=\tfrac12 S_t(Z_t)f_t\). It is pointwise continuous in \(t\), and (17), (19) yield
\[
|\partial_tf_t(z)|^2
\le3(1+|z|)^2e^{1+|z|}\phi(z).
\tag{20}
\]
This is integrable. Dominated convergence proves \(L^2\)-continuity of the derivative, DQM at each fixed point of the path, and the fundamental theorem of calculus in \(L^2\). The score is centered, also directly from \(\int-(\phi w_t)'=0\).

With the convention \(H(P,Q)^2=\int(\sqrt p-\sqrt q)^2\), Minkowski and (18) now give the finite, nonlinear path bound
\[
\boxed{H(P_{1,t},P_{1,0})
\le\tfrac12\int_0^{|t|}\|S_{\operatorname{sgn}(t)s}\|_2\,ds
\le20\epsilon_\delta|t|+10\delta t^2.}
\tag{21}
\]
No unspecified higher-order likelihood remainder appears in this inequality.

## 6. Nonvanishing interaction derivative

**PROVED HERE.** Differentiation under the expectations is justified by (7), the bounded coefficient path, and the finite moment \(m_h=Ee^{h(Z)}\). The latter follows from the saved family bound with \(2\delta\sqrt\theta<1\).
The exact expression is
\[
I_\delta'(t)=b_\delta e^{a_0+t}
\left\{E e^Z[
e^{a_2+u(Z+t+a_2)}-e^{u(Z+t)}]
-m_h(e^{b_2}-1)\right\}.
\tag{22}
\]
Oddness of \(h\), and \(h(z)\ge z\) for \(z\ge0\), imply
\(m_h=E\cosh(h(Z))\ge m_0=e^{1/2}\).
The bracket in (22) is at most
\[
m_0[e^{.7}-e^{-.1}-e^{.9}+1]\le-m_0/10.
\]
This sign and constant require no quadrature:
\(e^{.9}-e^{.7}=\int_{.7}^{.9}e^x\,dx\ge.2\) and
\(e^{-.1}\ge.9\).
Since \(b_\delta e^{a_0+t}\ge e^{-.1+.75-.05}=e^{.6}\),
\[
\boxed{I_\delta'(t)\le-k,\qquad
k=e^{.6}m_0/10=e^{11/10}/10.}
\tag{23}
\]
Thus \(|I_\delta(t)-I_\delta(s)|\ge k|t-s|\) throughout the path. The target does not flatten together with the observed score.

The [lower bound](lower_bound.md) applies (18), (21), (23) to the complete five-group experiment. All profile and coefficient assumptions have been checked before making that inference.

## 7. Explicit checks added in the v2 continuation

**PROVED HERE — verification details.** These calculations expand steps in the retained proof; they do not alter its constants or conclusion.

For the saturation tail bound, take \(v>A\) and set
\(a=(v+A)/(2\delta)\), \(b=(v-A)/(2\delta)>0\).
Then
\[
L'(v)=\tfrac12(\tanh a-\tanh b)
\le\tfrac12(1-\tanh b)\le e^{-(v-A)/\delta},
\]
\[
|L''(v)|=\frac{\operatorname{sech}^2 b-\operatorname{sech}^2 a}{4\delta}
\le\delta^{-1}e^{-(v-A)/\delta}.
\]
Evenness of \(L'\) and oddness of \(L''\) cover \(v<-A\). For \(|v|\le A\), the two squared hyperbolic secants lie in \([0,1]\), so their difference has magnitude at most one. This proves the global \(\delta^{-1}\chi\) bound without adding two tail bounds at the transition.

For \(r\ge A\) and \(j=1,2\),
\[
\frac{d}{dr}[r e^{-j(r-A)/\delta}]
=e^{-j(r-A)/\delta}(1-jr/\delta)\le0.
\]
Thus \(\sup_{r\ge0}r\chi(r)^j\le A\). In particular
\[
u'^2\le\chi(v)^2[\delta^2+4\delta\sqrt\theta|v|]
\le\delta^2+4\delta\sqrt\theta A,
\]
\[
|u''|\le\chi(v)[\delta+4\sqrt\theta|v|+2\delta\sqrt\theta]
\le\delta+4\sqrt\theta A+2\delta\sqrt\theta.
\]
The moving saturation zone therefore obeys the same original derivative envelopes as the central region, even though \(L''\) separately has a factor \(\delta^{-1}\).

Here is the central score arithmetic before rounding constants. Put \(q=|z|\) and \(E=E_\delta\). Equations (14)–(15) imply
\[
|w|\le3\delta|t|(1+q)+\tfrac32\delta E,
\]
\[
|S_t|\le\delta|t|[3q(1+q)+9/4+9\delta(1+q)]
+\delta E[(3/2)q+6(1+q)+(9/2)\delta].
\]
Since \(9\delta<1\), these bracketed polynomials are bounded by \(4(1+q)^2\) and \(8(1+q)\). Consequently the central norm is at most
\(20\delta|t|+16\delta E\), which is stronger than (16). The Gaussian tail bound gives exactly
\[
9e^{-r_\delta^2/8}=9e^{-A/(128\delta)}.
\]
Using the conservative (16), the total norm is at most
\(30\delta|t|+39\epsilon_\delta\), proving the stated constant 40. No unspecified polynomial absorption or uncomputed threshold in \(\delta\) is needed.

For (20), the flow bound implies
\(1+|Z_t(z)|\le(11/10)(1+|z|)\). Hence
\[
\tfrac14|S_t(Z_t)|^2\le\tfrac94(11/10)^2(1+|z|)^2
=\tfrac{1089}{400}(1+|z|)^2<3(1+|z|)^2.
\]
Together with (19), this gives the displayed integrable envelope. Also, changing variables \(x=Z_t(z)\) shows the exact identity
\[
\|\partial_t f_t\|_{L^2(dz)}^2
=\tfrac14\int S_t(x)^2\phi(x)\,dx.
\]
Thus the norm integrated in (21) really is the fixed-outcome density derivative in a common \(L^2\) space. The argument proves DQM on the interior of the full segment, with one-sided endpoint derivatives sufficient for its integral.

Finally, target differentiation can be justified here without importing an unverified moment claim. The warp envelope gives
\[
e^{h(z)}\phi(z)\le (2\pi)^{-1/2}
\exp\{-(1/2-\delta\sqrt\theta)z^2+|z|+\delta/(2\sqrt\theta)\}.
\]
The quadratic coefficient is strictly negative throughout the declared interval. Uniform bounded shifts and \(g_1'\le(5/4)e^x\) dominate the first-profile derivative by a constant times \(e^z\phi(z)\). These are integrable dominators for (22). The sign calculation is equivalently
\[
e^{.7}-e^{-.1}-e^{.9}+1
\le-\int_{.7}^{.9}1\,dx+(1-.9)=-.1.
\]
The target margin therefore survives the entire fixed-profile nonlinear path.
