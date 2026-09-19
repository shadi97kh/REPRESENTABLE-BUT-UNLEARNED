# Observable commutator, compact inversion and full-class calibration

All parameters and observations are as in [problem.md](problem.md). The algebra and constants below are analytic derivations; no numerical checks were run.

## 1. Exact six-term identity

**PROVED HERE.** For an increasing known onto warp let \(T=h^{-1}(h+1)\), \(T^{-1}=h^{-1}(h-1)\), and
\[
C(x)=T^{-1}(T(x+1)-1).
\]
The inverse is a function inverse, not a reciprocal. Set \(b_x=T(x+1)\). The observed anchor differences give
\[
R_1(x)=F_3(x)-F_0(x)=g_1(x+1)-g_1(x),
\quad U_1(x)=F_0(Tx)-F_4(x)=g_1(Tx)-g_1(x).
\]
Because \(T(Cx)=b_x-1\),
\[
R_1(x)+U_1(x+1)-R_1(b_x-1)-U_1(Cx)
=g_1(Cx)-g_1(x).
\]
Expanding the four observable expressions cancels \(F_0(b_x-1)\) and gives exactly
\[
\boxed{E(F;x)=F_3(x)-F_0(x)+F_0(b_x)-F_4(x+1)
-F_3(b_x-1)+F_4(Cx)=g_1(Cx)-g_1(x).}\tag{1}
\]
This uses only laws 0, 3 and 4. Replacing the profiles by any admissible profile derivatives \(r_1,r_2\) gives \(E(q;x)=r_1(Cx)-r_1(x)\). All four coefficient directions have zero anchor-law derivatives, so the same tangent identity holds for mixed paths. Evaluating an observed law at another quantile does not introduce another observed action.

## 2. A quantitative warp-class lemma

**PROVED HERE.** Suppose \(s\) is odd and \(C^3\), \(s'\ge0\) globally, \(|s'|,|s''|,|s'''|\le M\) on \([-3,3]\), \(M\ge1\), and \(s'''\ge\kappa>0\) on \([-2,2]\). Put \(h_\delta=\mathrm{id}+\delta s\). It is increasing and onto, with \(h_\delta'\ge1\). For
\[
0<\delta\le\min\{1/M,\kappa/(500M^2)\},\quad J=[-2,0],
\]
the composition in (1) obeys
\[
C_\delta(J)\subset J,\qquad
0<C_\delta'\le\lambda_\delta:=1-\kappa\delta/2<1. \tag{2}
\]
The proof verifies the proposed Taylor constants, rather than assuming them.

Let \(y=T_{\delta,t}(x)=h_\delta^{-1}(h_\delta(x)+t)\), \(t=\pm1\). Write
\[
B=1+\delta s'(y),\quad A=1+\delta s'(x),\quad
a=y_x,\ d=y_\delta,\ e=y_{\delta\delta}.
\]
Since \(h_\delta'\ge1\), \(|y-x|\le1\), with the sign of \(t\). On the coordinate domain below, \(1\le A,B\le2\) whenever \(\delta M\le1\). Implicit differentiation yields
\[
Ba=A,\quad Bd=s(x)-s(y),\quad
Be=-2s'(y)d-\delta s''(y)d^2,
\]
\[
By_{x\delta}=s'(x)-s'(y)a-\delta s''(y)da,\quad
By_{xx}=\delta[s''(x)-s''(y)a^2],
\]
\[
By_{xxx}=\delta[s'''(x)-s'''(y)a^3-3s''(y)ay_{xx}].
\]
Consequently
\[
|d|\le M,\quad |e|\le3M^2,\quad \tfrac12\le a\le2,\quad
|y_{x\delta}|\le5M,\quad |y_{xx}|\le5\delta M,\quad
|y_{xxx}|\le39\delta M.
\]
For completeness, the two remaining mixed derivatives follow from
\[
Be_x=-\delta s''(y)ae-2s''(y)ad-2s'(y)d_x
-\delta s'''(y)ad^2-2\delta s''(y)dd_x
\]
and
\[
By_{xx\delta}
=s''(x)-s''(y)a^2-\delta s'''(y)da^2
-2\delta s''(y)ay_{x\delta}
-[s'(y)+\delta s''(y)d]y_{xx}.
\]
Their termwise bounds are respectively
\[
|y_{x\delta\delta}|\le(6+4+10+2+10)M^2=32M^2,
\quad |y_{xx\delta}|\le(1+4+4+20+5+5)M=39M.
\]
No fourth derivative of \(s\) is required.

The coordinate domain is valid before any contraction conclusion. If \(x\in J\), then \(x+1\in[-1,1]\),
\[
T_{\delta,+1}(x+1)\in[-1,2],\quad
u:=T_{\delta,+1}(x+1)-1\in[-2,1],\quad
T_{\delta,-1}(u)\in[-3,1].
\]
The intervals between input and output of either flow are within \([-3,3]\). These inclusions hold for every intermediate parameter \(0\le\delta'\le\delta\).

Let \(f=T_{\delta,-1}\). Differentiating \(C_\delta(x)=f(u)\) twice in \(\delta\) gives
\[
\partial_\delta^2C=f_{\delta\delta}+2f_{x\delta}u_\delta
+f_{xx}u_\delta^2+f_xu_{\delta\delta}.
\]
The preceding inequalities bound this by \((3+10+5+6)M^2=24M^2\). Its \(x\)-derivative consists of
\[
f_{x\delta\delta}u_x+2f_{xx\delta}u_xu_\delta
+2f_{x\delta}u_{x\delta}+f_{xxx}u_xu_\delta^2
+2f_{xx}u_\delta u_{x\delta}
+f_{xx}u_xu_{\delta\delta}+f_xu_{x\delta\delta},
\]
bounded by \((64+156+50+78+50+30+64)M^2=492M^2\).
At \(\delta=0\), \(C_0(x)=x\) and
\[
\partial_\delta C_0(x)=2s(x+1)-s(x+2)-s(x)=-D_s(x).
\]
Integral Taylor remainder therefore proves
\[
C_\delta(x)=x-\delta D_s(x)+r_\delta(x),\qquad
\|r_\delta\|_\infty\le12M^2\delta^2,\quad
\|r_\delta'\|_\infty\le246M^2\delta^2. \tag{3}
\]

Oddness gives \(D_s(-1)=0\), while
\[
D_s'(x)=\int_0^1\!\int_0^1s'''(x+u+v)\,du\,dv\ge\kappa.
\]
Thus \(D_s(-2)\le-\kappa\), \(D_s(0)\ge\kappa\), and
\[
C_\delta(-2)>-2,\quad C_\delta(0)<0
\]
because \(12M^2\delta^2\le(12/500)\kappa\delta\).
Both flows are increasing. Finally,
\[
C_\delta'\le1-\kappa\delta+246M^2\delta^2
\le1-\kappa\delta/2.
\]
This proves (2). The same flow derivative ratios give \(C_\delta'\ge1/4\); only its positivity is needed statistically. In particular \(r=-1\) need not be a fixed point of \(C_\delta\).

## 3. Verification for the explicit family

**PROVED HERE.** Direct differentiation gives
\[
s_\theta'=\frac{1+2\theta x^2}{\sqrt{1+\theta x^2}},\quad
s_\theta''=\frac{\theta x(3+2\theta x^2)}{(1+\theta x^2)^{3/2}},
\quad s_\theta'''=\frac{3\theta}{(1+\theta x^2)^{5/2}}.
\]
For \(\theta\in[1/2,2]\) and \(|x|\le3\),
\[
|s_\theta'|\le2\sqrt{19}<10,\quad
|s_\theta''|\le3\sqrt2<10,\quad |s_\theta'''|\le6<10.
\]
The second bound follows by writing \(u=\sqrt\theta |x|\) and bounding
\(u(3+2u^2)/(1+u^2)^{3/2}\le3\).
On \([-2,2]\), the minimum of \(3\theta/(1+4\theta)^{5/2}\) occurs at \(\theta=2\), since its derivative has sign \(1-6\theta\). Hence
\[
s_\theta'''\ge2/81>1/64.
\]
The public choices \(M=10,\kappa=1/64\) therefore imply precisely
\[
\boxed{0<\delta\le1/3{,}200{,}000.}\tag{4}
\]
This proof supplies no conclusion for \(\delta=.1\) from (4).

After invariance, all six source-query arguments in (1) belong to \([-2,2]\): \(x,x+1,b_x,b_x-1,Cx\). Each query map has Lipschitz constant at most 2 on \(J\). These are only the **central** queries; calibration below uses much larger arguments.

## 4. Finite normalized reconstruction and transports

**PROVED HERE.** Put \(r=-1\), and
\[
D_K(F;x,r)=-\sum_{k=0}^{K-1}[E(F;C^kx)-E(F;C^kr)].
\]
Telescoping (1) proves the exact finite identity
\[
g_1(x)=1-[F_3(r)-F_0(r)]+D_K(F;x,r)
+g_1(C^Kx)-g_1(C^Kr). \tag{5}
\]
For \(x\in J\), the final difference has absolute value at most
\(2\lambda_\delta^K|x-r|\le2\lambda_\delta^K\), since \(g_1'\le2\) on \(J\).
The two iterates must be retained; replacing either with an assumed known fixed-point profile value is invalid.

For any specified integer chart \(m\) with \(x+m\in J\), define
\[
H_m(F;x)=
\begin{cases}
-\sum_{j=0}^{m-1}(F_3-F_0)(x+j),&m>0,\\
\sum_{j=1}^{-m}(F_3-F_0)(x-j),&m<0,\\
0,&m=0.
\end{cases}
\]
The finite reconstruction is
\[
G_{1,K}(F;x;m)=1-[F_3(r)-F_0(r)]
+D_K(F;x+m,r)+H_m(F;x).
\]
At population curves its error is just the terminal difference in (5) at \(x+m\); no error proportional to transport length is added to the population identity. Set
\[
G_{2,K}(F;y;m)=F_0(h^{-1}y)-G_{1,K}(F;h^{-1}y;m).
\]
For arbitrary isolated evaluations a public integer may be chosen to move the argument into \([-1,0]\). For moving calibration evaluations use the fixed whole-interval charts below; an empirical reconstruction need not be globally continuous across changing charts.

The [conditioning proof](conditioning.md) proves a kernel increment bound \(C_\kappa|x-y|^{1/2}\) and therefore the Hilbert-space series bound
\[
\|d_{x,r}\|_\pi\le
\frac{6\sqrt2 C_\kappa}{1-\sqrt{\lambda_\delta}}\sqrt{|x-r|}
\le\frac{24\sqrt2 C_\kappa}{\kappa\delta}\sqrt{|x-r|}. \tag{6}
\]
This establishes bounded central profile-evaluation functionals, after adding the normalizer, and bounded finite transported evaluations. It is not the complete target variance or the efficient bound.

## 5. Coefficients over the original profile class

**PROVED HERE.** Choose the positive known \(R\) solving \(\delta s_\theta(R)=3\). Squaring \(R\sqrt{1+\theta R^2}=3/\delta\) yields
\[
R^2=\frac{\sqrt{1+36\theta/\delta^2}-1}{2\theta},
\qquad h(R)=R+3,\quad h(-R)=-R-3. \tag{7}
\]
For either coefficient pair, with \(w_-=e^R,\ w_+=e^{-h(R)}\), set
\[
B(c)=
\begin{pmatrix}
w_-[g_1(-R+\alpha)+g_2(-h(R)+\beta)]\\
w_+[g_1(R+\alpha)+g_2(h(R)+\beta)]
\end{pmatrix}.
\]
The observed counterpart is \(y_i=(w_-F_i(-R),w_+F_i(R))^\top\).
On the original box, the Jacobian diagonal entries are at least \(\frac12e^{.5}\); each off-diagonal entry is at most \(2e^{1-3}\). This follows by applying the original derivative envelope separately to all four entries, with the indicated row scalings.

For the enlarged search box \(\mathcal R=[.4,1.1]^2\), define
\[
a_0=\tfrac12e^{.4},\quad A_0=2e^{1.1},\quad b_0=2e^{-1.9},
\quad\gamma=a_0-b_0>0,\quad \tau=A_0^{-1},\quad q=1-\tau\gamma<1. \tag{8}
\]
Positivity of \(\gamma\) is equivalent to \(e^{2.3}>4\), already implied by the first three positive terms of the exponential series.
Every diagonal entry lies in \([a_0,A_0]\), and every off-diagonal lies in \((0,b_0]\). Averaging the Jacobian along any coefficient segment preserves these inequalities. In a row corresponding to the largest component of \(c-c'\), the reverse triangle inequality gives
\[
\|B(c)-B(c')\|_\infty\ge\gamma\|c-c'\|_\infty,
\qquad \|J(c)^{-1}\|_\infty\le\gamma^{-1}. \tag{9}
\]
This proves global separation, not just local Jacobian invertibility.

Since \(1-\tau J_{jj}\ge0\), the absolute row sums of \(I-\tau J(c)\) are at most \(1-\tau a_0+\tau b_0=q\). Segment integration and nonexpansiveness of coordinatewise projection prove that
\[
c\longmapsto \Pi_{\mathcal R}\{c+\tau[y_i-B(c)]\} \tag{10}
\]
is a population contraction in infinity norm. Its fixed point is the true pair. The step size uses no true-profile Jacobian. The profile envelope has not been narrowed by a factor depending on \(\delta\).

## 6. Fixed charts for every calibration query

**PROVED HERE.** The four public latent-coordinate query ranges are
\[
\mathcal I_{1,\pm}=\pm R+[.4,1.1],\qquad
\mathcal I_{2,\pm}=h^{-1}(\pm h(R)+[.4,1.1]).
\]
Each is a closed interval \([L_j,U_j]\) of width at most \(7/10\), since \(h^{-1}\) is 1-Lipschitz. Choose once, from known endpoints,
\[
m_j=\lceil-7/4-L_j\rceil.
\]
Then
\[
[L_j+m_j,U_j+m_j]\subset[-7/4,-1/20]\subset(-2,0). \tag{11}
\]
Indeed \(L_j+m_j\in[-7/4,-3/4)\), and adding at most \(7/10\) gives an upper endpoint strictly below \(-1/20\). Both profiles and both signs therefore have a fixed chart on the entire search box, with positive margins. This uses the width-two overlap of \(J\).

All original latent arguments have absolute value at most \(R+1.1\). Their finite transports and normalization lie in \([-D,D]\), \(D=R+3\), and \(|m_j|\le L:=\lceil R+4\rceil\). These crude bounds follow because a unit transport stays between its original argument and its endpoint in \(J\). The central commutator queries remain in \([-2,2]\).

An exact ceiling can be avoided in certified arithmetic: enclose each known \(L_j\) within width \(1/100\), choose \(m_j\) by its rational lower endpoint, and then the transported upper endpoint remains below \(-1/25\). The proofs use only fixed charts with these margins, not the particular ceiling convention. This is an explicitly permitted chart-selection variant with its corresponding influence, not an assumption that exact transcendental integer comparisons can always be decided.

**OPEN:** diagonal dominance alone says nothing favorable about estimating \(F_i(\pm R)\). The radius grows like \(\delta^{-1/2}\), and the individual quantile influences incur Gaussian-tail costs. They are explicitly retained in [conditioning.md](conditioning.md). The [attainment proof](attainment.md) uses population contraction plus empirical equicontinuity; it never assumes the empirical update contracts.
