> Publication copy of `docs/interaction_recoverability/proofs.md`. Mathematical and result content is retained; execution material is omitted. This dated record retains its original conclusions and limitations, including conclusions superseded by later stages. Original SHA-256: `350fb7517f14efe6aab31d724bdf0f8d15b158688dd24d42b2b45d6a6b7afd0f`.

# Proof notes and boundaries of interaction recoverability

Each claim below has an explicit status. Numerical checks support implementation only. The rate derivation uses established likelihood, two-point and empirical-CDF tools and is classified as incremental; these notes do not claim a new general theorem on unknown latent mechanisms.

**PROVED HERE — validity and integration.** Put \(s(z)=z\sqrt{1+z^2}\). Then

\[
s'(z)=\frac{1+2z^2}{\sqrt{1+z^2}}>0,\quad
s''(z)=\frac{z(3+2z^2)}{(1+z^2)^{3/2}},\quad |s''(z)|\le2.
\]

The last bound follows, for example, from \(s'''(z)=3/(1+z^2)^{5/2}>0\) and limits \(s''(\pm\infty)=\pm2\). Thus \(h_\delta'=1+\delta s'>0\), \(h_\delta(0)=0\), and \(h_\delta\) is real analytic and onto the real line. Positive sums \(A e^z+B e^{h_\delta(z)}\) increase strictly from zero to infinity. Both laws consequently have the same support \((0,\infty)\), with positive densities and no endpoint atoms.

For \(0\le\delta\le.1\), \(|s(z)|\le z^2+1/2\) gives the integrable envelope

\[
e^{h_\delta(z)}\phi(z)\le(2\pi)^{-1/2}
\exp\{-0.4z^2+|z|+0.05\}.
\]

Bounded action coefficients multiply this by only a fixed constant. Dominated convergence is therefore justified for response means and absolute differences. More generally \(E e^{p h_\delta(Z)}<\infty\) for \(p\delta<1/2\); at equality and above the positive tail diverges. Variances exist throughout the specified class, but at \(\delta=.1\) the fifth moment does not. Neither bounded outcomes nor arbitrary exponential-moment assumptions are available.

**PROVED HERE — axes, quantile coupling and limits.** At \(\delta=0\), T and Q agree pointwise on each axis. At positive delta they still agree pointwise on axes 1, 3 and 4. At \(t e_2\),

\[
T-Q=(e^t-e^{t/2})(e^{h_\delta(z)}-e^z).
\]

For \(t\in[0,1]\), its W1 is the absolute coefficient times \(E|e^{h_\delta(Z)}-e^Z|\), tending to zero. This describes all points on the second axis; the statistical experiment below observes only its specified discrete stratum. Because both transformations increase strictly in the same continuously distributed scalar latent variable, their shared latent coupling is the common-quantile coupling and is optimal for scalar W1. This assertion would generally fail for nonmonotone or multivariate generators.

At \(a_{12}=(1,1,0,0)\), the zero-delta pointwise difference is \((e+e^2-2e^{3/2})e^z>0\). Consequently its limiting W1 is

\[
L=e^{1/2}(e+e^2-2e^{3/2})=1.886070833180237\ldots.
\]

The primary *interaction* gap is not the positive-delta joint W1. It is exactly

\[
I_{12}(T)-I_{12}(Q)
=(e-\sqrt e)\{(e-1)m_\delta-(\sqrt e-1)m_0\}\ge L.
\]

To prove the inequality, pair positive and negative z in \(m_\delta\). Their sum is proportional to \(2\cosh(z+\delta s(z))\), which increases with delta for positive z; hence \(m_\delta\ge m_0\). This gives a nonvanishing target separation for every allowed delta, rather than just its limit.

**STANDARD RESULT APPLIED — relation to DExtrI.** In each fixed positive-delta two-model class, zero additive functions and exponential profiles satisfy the profile/analytic requirements; the two distinct supports are nontrivial; the median-zero increasing warps are unbounded and nonproportional; and \(h_\delta(z)/z\to\infty\). These check the relevant A1–A7 structural conditions. Scale labels are fixed as in the problem specification. The theorem's interval observation design is a separate requirement and is absent in our finite dictionary. At zero delta, nonproportionality fails. A union across positive deltas also fails the pooled distinct-growth condition: \(h_\delta/h_{\delta'}\to\delta/\delta'\), neither zero nor infinity. Thus this is not a counterexample to a population identification theorem. [DExtrI v1, Assumptions 3.2 and Theorems 3.3–3.4](https://arxiv.org/html/2608.19849v1).

**PROVED HERE — density and uniform transport-score bound.** Consider \(F_t(z)=A e^z+B e^{z+t s(z)}\), with \(A,B>0\) fixed and \(t\in[0,.1]\). Its density is \(p_t(y)=\phi(z)/F'_t(z)\), where \(z=F_t^{-1}(y)\). Define \(u_t=\partial_tF_t/F'_t\). Differentiating at *fixed y*, including the inverse map and Jacobian, gives

\[
\partial_t\log p_t(F_t(z))=z u_t(z)-u_t'(z).
\]

With \(r=(B/A)e^{ts}\),

\[
u_t=\frac{s}{r^{-1}+1+t s'},\qquad
|u_t|\le|s|,\quad |u_t'|\le s'+t|s|(s'+2).
\]

For the derivative bound, differentiate the denominator: its derivative is \(-t s' r^{-1}+t s''\). Divide by the denominator, use its lower bound one and \(r^{-1}/(r^{-1}+1+t s')\le1\). Since \(s'\le1+2|z|\), the score is bounded by

\[
P(|z|)=1.2|z|^3+0.3|z|^2+2.6|z|+1.15.
\]

Normal absolute moments give the explicit uniform information bound

\[
\int (\partial_t\log p_t)^2p_t\le E P(|Z|)^2
=C_0=65.6233873491624\ldots.
\]

This bound is independent of A and B. On the common support, inverse functions are continuously differentiable in t, so pointwise differentiation of \(\sqrt{p_t}\) and the fundamental theorem of calculus apply. Its L2 derivative norm is at most \(\sqrt{C_0}/2\); Minkowski's integral inequality then gives \(H(P_t,P_0)\le t\sqrt{C_0}/2\). This argument controls tails through the normal score moment and does not assume a numerical density integration is a proof.

**STANDARD RESULT APPLIED — complete-experiment lower bound for the supplied pair.** Use the convention \(H^2(P,Q)=\int(\sqrt p-\sqrt q)^2=2-2\rho(P,Q)\). T and Q at \(e_2\) share their zero-delta law. The preceding two paths and the triangle inequality give \(H^2(P_{T,e_2},P_{Q,e_2})\le C_0\delta^2\). All other observed strata have identical laws. Therefore

\[
H^2(\mathbb P_T,\mathbb P_Q)
=2\{1-(1-H^2(P_{T,e_2},P_{Q,e_2})/2)^{n_{e_2}}\}
\le C_0 n_{e_2}\delta^2.
\]

Since \(\operatorname{TV}\le H\), thresholding any contrast estimator halfway between the two target values and using the testing error bound yields

\[
\inf_{\widehat I}\max_{M=T,Q}E_M(\widehat I-I_M)^2
\ge \frac{L^2}{8}\big(1-\min\{1,\delta\sqrt{C_0 n_{e_2}}\}\big).
\]

For \(\delta\sqrt{C_0 n_{e_2}}\le1/2\), this is at least \(L^2/16\). Counts at the four identical strata cannot change this pairwise difficulty. This is a standard two-point argument applied to the *entire* observed product law, not an implication from Wasserstein proximity.

**PROVED HERE — amplitude path for the larger normalized class.** Let \(m=\sqrt e\), \(D=e-m\), and hold the coefficients at action \(e_1\) equal to \((m,e)\). At \(e_2\), take \((u_2,v_2)=(b,m+e-b)\). Vary b within \([m,e]\), leaving every other observed stratum unchanged. For \(F_b(z)=b e^z+(m+e-b)e^{h_\delta(z)}\), its transport velocity is

\[
w_b(z)=\frac{1-e^{\delta s(z)}}{b+(m+e-b)e^{\delta s(z)}(1+\delta s'(z))}.
\]

Using \(|1-e^t|\le |t|\max(1,e^t)\), both coefficients at least m, and the same denominator-derivative argument gives

\[
|w_b|\le\delta|s|/m,\qquad
|w_b'|\le(\delta/m)\{s'+.1|s|(s'+2)\}.
\]

Thus its b-score has squared mean at most \(\delta^2 C_0/m^2\), uniformly in the whole path. Two centered values separated by \(\eta\le D\) have single-observation Hellinger distance at most \(\eta\delta\sqrt{C_0}/(2m)\). The target derivative is \(m_0(m-1)-m_\delta(e-1)\), whose magnitude is at least \(\kappa=m_0D\). Choose

\[
\eta=\min\{D,m/(\delta\sqrt{C_0 n})\}.
\]

Their complete balanced experiment has TV at most 1/2 and target separation at least \(\kappa\eta\).

**STANDARD RESULT APPLIED — minimax lower rate.** The same testing reduction proves

\[
R_n(I_{12},\mathcal M_\delta)\ge
\frac{\kappa^2}{16}\min\{D^2,m^2/(C_0n\delta^2)\}.
\]

No local asymptotic approximation or assumption of Fisher nonsingularity is used. This lower bound also covers procedures that use all five strata. It extends the fixed T/Q example to the high-information regime, where a fixed two-point separation alone would give a vacuous bound.

**STANDARD RESULT APPLIED — constructive upper rate.** At each of \(e_1,e_2\), the .25 and .75 response quantiles equal \(V_\delta(u_i,v_i)^\top\), where \(r=\Phi^{-1}(.75)\),

\[
V_\delta=\begin{pmatrix}e^{-r}&e^{h_\delta(-r)}\\e^r&e^{h_\delta(r)}\end{pmatrix},
\quad \det V_\delta=2\sinh(\delta s(r)).
\]

Consequently \(\|V_\delta^{-1}\|_\infty\le C_V/\delta\) uniformly over allowed delta, from bounded matrix entries and \(\sinh t\ge t\). This is an observable-central-quantile separation under a known basis, not an assumed extrapolation inequality.

Let \(D_n\) be the maximum empirical-CDF sup error over these two strata. The [DKW–Massart inequality](https://arxiv.org/html/2607.04387v1) and a union bound give \(P(D_n>t)\le4e^{-2nt^2}\). For \(D_n\le1/8\), the four empirical quantile levels lie in [.125,.875]. Set \(r_*=\Phi^{-1}(.875)\). The known uniform quantile Lipschitz bound is

\[
L_\delta=\frac{e\{e^{r_*}+e^{h_\delta(r_*)}(1+\delta s'(r_*))\}}{\phi(r_*)}.
\]

It is bounded uniformly in delta. Invert the empirical quantile vectors, clip coefficients to \([m,e]\), and substitute in the explicit formula for \(I_{12}\). Clipping cannot increase coordinate error. With \(K_\delta=2(e-1)(m_0+m_\delta)\), target error is at most \(K_\delta\) times the maximum coefficient error. Integrating the DKW tail gives \(E D_n^2\le(\log4+1)/(2n)\). Hence a finite bound is

\[
E(\widehat I-I)^2\le \min\left\{K_\delta^2D^2,
\frac{K_\delta^2L_\delta^2\|V_\delta^{-1}\|_\infty^2(\log4+1)}{2n}
+4K_\delta^2D^2e^{-n/32}\right\}.
\]

All constants except the displayed inverse factor remain uniformly bounded. Since \(\sup_{n\ge1}ne^{-n/32}<\infty\) and \(\delta^2\le.01\), the exponential remainder is absorbed into a constant times \((n\delta^2)^{-1}\). This proves the matching rate order stated in the problem file. It does not establish good constants, likelihood efficiency, or a novel estimator. Known quadrature error in \(m_\delta\) would add a bounded deterministic plug-in error; the mathematical theorem defines that moment exactly.

**PROVED HERE — a target that avoids the weak direction.** Pointwise,
\(G(e_1+e_3+e_4,z)=eG(e_1,z)\) and \(G(e_3+e_4,z)=eG(0,z)\). Therefore the secondary J equals \((e-1)(\mu(e_1)-\mu(0))\). Its sample-mean estimator is unbiased with variance \((e-1)^2\{\operatorname{Var}(Y_{e_1})/n_{e_1}+\operatorname{Var}(Y_0)/n_0\}\). The variance bound is uniform in delta by the second-moment calculation above. This is ordinary functional estimability under a known operator; it is not a new general principle.

**PROVED HERE — unrestricted misspecification destroys the conclusion.** Adding \(\eta a_1a_2\) to the response leaves every observed law unchanged and shifts \(I_{12}\) by eta. If \(\eta\in[-B,B]\), two indistinguishable endpoints separated by 2B imply minimax MSE at least \(B^2\): the average of the two squared errors is at least \(B^2\). If eta is unrestricted, finite uniform risk is impossible. This perturbation lies outside the normalized exponential class; calling it a robustness guarantee would be false.

**OPEN — exact remaining proof gap.** No result here estimates unknown profiles/warps or a separation constant from the same limited action dictionary while retaining the displayed target rate. Nor is there a result that adapts to unknown misspecification magnitude or supplies sharper target-specific inference than existing inverse-functional theory. A credible extension must formulate the nuisance class and measurable central-range design condition, establish a target modulus without assuming the answer, and construct an estimator attaining it with nuisance-estimation error included. The restricted known-warp lower and upper proofs above have no intentionally unproved lemma; extending them to that broader claim is the open step. The novelty of such an extension is also unresolved until its precise statement is compared with the sources.
