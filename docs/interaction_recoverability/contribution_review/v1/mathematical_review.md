# Combined mathematical review — 13 September 2026

This internal static assessment supports a model-specific constructive attainment theorem paired with an exponential worst-case efficient-variance obstruction and an explicit finite-experiment risk lower bound. The corollary below rules out extending the historical exponential-profile risk order uniformly to the original profile class. No fatal defect was identified in the retained proof chain during this review; that is an internal assessment, not independent mathematical certification or external review. The pointwise upper theorem remains the principal technical dependency.

The [dependency map](dependency_map.md) identifies exact, unchanged proof copies and separates mathematical dependencies from historical and implementation evidence. [The decision](decision.md) records the claims left open. All calculations in this document are analytic; no scientific program was run.

## 1. One experiment and three different conclusions

Fix **known** \(\theta\in[1/2,2]\) and **known** \(0<\delta\le\delta_*:=1/3{,}200{,}000\). Set
\[
h(z)=z+\delta s_\theta(z),\qquad s_\theta(z)=z\sqrt{1+\theta z^2}.
\]
Write \(\eta=(g_1,g_2,\alpha_1,\beta_1,\alpha_2,\beta_2)\in\Theta\), where all four coefficients are unknown in \([.5,1]\), and each profile belongs to the **full original** envelope
\[
\mathcal G=\{g:\mathbb R\to(0,\infty):g\text{ real analytic},\ g(0)=1,\ g(-\infty)=0,
\quad \tfrac12 e^x\le g'(x)\le2e^x,\quad |g''(x)|,|g'''(x)|\le10e^x\}.
\]
Real analytic does not mean complex entire. The profile class is not a neural family or a neighborhood of the exponential pair.

The data are five independent lists \(Y_{ai}=F_a(Z_{ai})\), \(a=0,\ldots,4\), \(i=1,\ldots,n\), with independent, unobserved \(Z_{ai}\sim N(0,1)\). Thus **\(N=5n\)**. Within one response both summands share its latent draw; observations across actions have no latent pairing. The five increasing curves are
\[
\begin{aligned}
F_0(z)&=g_1(z)+g_2(h(z)),\\
F_i(z)&=g_1(z+\alpha_i)+g_2(h(z)+\beta_i),&&i=1,2,\\
F_3(z)&=g_1(z+1)+g_2(h(z)),\\
F_4(z)&=g_1(z)+g_2(h(z)+1).
\end{aligned}
\]
The component labels and unit private shifts in sources 3 and 4 are known. The unavailable joint-action law has shifts \((\alpha_1+\alpha_2,\beta_1+\beta_2)\). For \(\Delta_{u,v}g(x)=g(x+u+v)-g(x+u)-g(x+v)+g(x)\), the target and risk are
\[
I_\eta=E\Delta_{\alpha_1,\alpha_2}g_1(Z)+E\Delta_{\beta_1,\beta_2}g_2(h(Z)),\qquad
\mathcal R_{N,\delta,\theta}=\inf_{\widehat I(Y_0,\ldots,Y_4)}\sup_{\eta\in\Theta}E_\eta(\widehat I-I_\eta)^2.
\]
All expectations here are on the declared response scale. The infimum may include randomized estimators; the testing proof still applies.

The combined statement has three parts:

1. **Pointwise regular attainment.** The specified empirical-quantile statistic satisfies, at every fixed \((\delta,\theta,\eta)\),
   \[
   \widehat I_N-I_\eta=\frac1N\sum_{a=0}^4\sum_{i=1}^n\psi_{\eta,a}(Y_{ai})+o_{P_\eta}(N^{-1/2}),\qquad E_a\psi_{\eta,a}=0.
   \]
   Its complete influence is square integrable and gives a root-\(N\) CLT, regular along each fixed admissible DQM path. Its variance need not be efficient. This is a mathematical statistic on exact observations; the optional certified arithmetic realization requires refinable observation digits and certified known primitives. It is not an implementation or a theorem for arbitrary fixed-resolution assay measurements.
2. **Worst-case efficient variance.** Let \(\|f\|_\pi^2=\frac15\sum_a E_af_a^2\). Define \(V_{\rm eff}\) using the minimum-norm compatible influence in the closure of actual admissible scores. Then
   \[
   \frac{e^{11/5}}{32000}e^{1/(640\delta)}\le\sup_{\eta\in\Theta}V_{\rm eff}(\eta;\delta,\theta)\le\Psi(\delta,\theta)^2,
   \]
   where \(\Psi\) is the public, truth-uniform constructed-influence bound in [conditioning, (4)–(10)](evidence/docs/interaction_recoverability/warp_family/v1/conditioning.md). Consequently
   \[
   \log\!\left(1+\sup_{\eta\in\Theta}V_{\rm eff}(\eta;\delta,\theta)\right)=\Theta(1/\delta)
   \]
   as \(\delta\downarrow0\), with constants uniform in \(\theta\in[1/2,2]\). This matches logarithmic order, not leading exponential constants.
3. **Finite-experiment squared-risk lower bound.** For every integer \(n\ge1\), with \(N=5n\),
   \[
   \boxed{\mathcal R_{N,\delta,\theta}\ge C_*\min\left\{1,\frac{e^{1/(640\delta)}}N,\frac1{\delta\sqrt N}\right\},\qquad C_*:=\frac{e^{11/5}}{2{,}048{,}000}.}
   \]
   This inequality has no asymptotic remainder and applies to shrinking-separation sequences in its declared range. There is no matching finite-sample MSE upper theorem.

The tangent convention matters: profile directions satisfy \(r_j(0)=0\), \(|r_j^{(k)}(x)|\le C_v e^x\), \(k=0,1,2,3\), and paths must remain admissible. At boundaries only genuinely admitted paths are asserted. The retained direction norm is
\[
\|v\|_X^2=|d_c|^2+\sum_{j=1}^2\sum_{k=0}^3\int |r_j^{(k)}(x)|^2\frac{e^{-2x}}{1+x^2}\,dx.
\]
The score closure is not replaced by any computed finite dictionary. For a curve velocity \(q_a=\dot F_a\), its actual density score in latent coordinates is \(S_a=zq_a/F_a'-(q_a/F_a')'\).

## 2. What the constructive upper proof actually supplies

The [commutator proof](evidence/docs/interaction_recoverability/warp_family/v1/commutator_and_calibration.md) defines \(T=h^{-1}\circ(h+1)\), \(C=T^{-1}(T(x+1)-1)\), \(b_x=T(x+1)\), and the observable identity
\[
F_3(x)-F_0(x)+F_0(b_x)-F_4(x+1)-F_3(b_x-1)+F_4(Cx)=g_1(Cx)-g_1(x).
\]
Cancellation can be checked directly from \(h(Tx)=h(x)+1\). On \(J=[-2,0]\), the proved warp bounds give \(0<C'\le\lambda=1-\delta/128<1\). Subtracting the telescoping series at \(x\) and \(r=-1\) reconstructs profile differences. The known normalization gives \(g_1(-1)=1-(F_3-F_0)(-1)\). The anchor point \(-1\) is not assumed to be a fixed point of \(C\). Finite unit transports recover other required evaluations; fixed charts prevent integer transport changes at moving calibration arguments from being silently differentiated.

Coefficients are calibrated at \(\pm R\), where \(\delta s_\theta(R)=3\) and
\[
R^2=\frac{\sqrt{1+36\theta/\delta^2}-1}{2\theta}.
\]
The scaled two-row calibration map has inverse Jacobian norm at most \(\gamma^{-1}\), \(\gamma=\tfrac12e^{.4}-2e^{-1.9}>0\), on the **computational** box \([.4,1.1]^2\). This box leaves all original truths interior to the search region. The contraction is a population comparison; it does not require the piecewise-linear empirical map itself to have a contracting derivative.

The [attainment specification](evidence/docs/interaction_recoverability/warp_family/v1/attainment.md), Sections 1–6, sorts each observed list and interpolates through \((i/n,Y_{a(i)})\); its quantile curves are evaluated at \(\Phi(z)\). It uses the explicit schedules
\[
K_N=\lceil2\log N/|\log\lambda|\rceil+2,\quad B_N=\sqrt{(3/2)\log N},\quad M_N=\lceil4B_N+30\rceil+2,
\]
\[
k_N=\lceil\log N/|\log q|\rceil,\qquad q=1-\gamma/(2e^{1.1}),
\]
with \(2k_N\) projected coefficient updates. The signed target identity retains periodized kernels, outward transport weights and their joins. The statistic includes finite tail/series truncations, a public exceptional-data branch and a bounded final clip. It never receives true profiles, latent values, calibration Jacobians, influences or joint outcomes.

The proof must control accumulated compact quantile remainders, their moving-query modulus, coefficient localization before linearization, growing central series, target tails, periodization and shared-source covariance. These are supplied in Sections 3–6, rather than inferred from a fixed-quantile CLT alone. Sections 7–8 state the certified arithmetic assumptions and tolerances. They do not license a floating-point backend.

The complete influence includes estimated coefficients and all reference-law reuse. In the [conditioning proof](evidence/docs/interaction_recoverability/warp_family/v1/conditioning.md),
\[
\|\psi_\eta\|_\pi\le\Psi=D_t+8D_c+B_I Z_R/\gamma,
\]
with every quantity explicit in (4)–(10). In particular \(D_c=O(\delta^{-1})\), whereas \(Z_R\) contains \(e^R Q(R+3)\) and \(Q(z)\) contains Gaussian extreme-quantile amplification \(e^{z^2/4}\). It follows that
\[
\log(1+\Psi^2)\le\frac3{2\delta\sqrt\theta}+O(\delta^{-1/2})+O(\log(1/\delta)),
\]
uniformly in \(\theta\). Compact central inversion does not remove remote calibration cost. Square-integrability and the pointwise remainder bound do not imply uniform integrability of \(N(\widehat I-I)^2\), nor an MSE upper bound. Even a bounded statistic and a pointwise CLT are insufficient for that inference.

## 3. Check of the fixed-profile hard path

The exact definitions are in [construction](evidence/docs/interaction_recoverability/warp_family_lower_bound/v2/construction.md); all global estimates and DQM details are in [admissibility and score](evidence/docs/interaction_recoverability/warp_family_lower_bound/v2/admissibility_and_score.md). Put
\[
A=1/10,\quad a_0=3/4,\quad a_2=3/5,\quad b_2=9/10,\quad |t|\le1/20,
\]
\[
L_\delta(v)=\delta\left[\log\cosh\frac{v+A}{2\delta}-\log\cosh\frac{v-A}{2\delta}\right],\quad
v(z)=\delta s_\theta(z),\quad u(z)=L_\delta(v(z)),
\]
\[
b=\left[\int_{-\infty}^0 e^{x+u(x-a_0)}dx\right]^{-1},\qquad
g_1(x)=b\int_{-\infty}^x e^{y+u(y-a_0)}dy,\quad g_2(x)=e^x.
\]
For each known \((\delta,\theta)\), freeze **these same profiles**, including their internal \(a_0\), throughout the path. Change only
\[
\alpha_1(t)=a_0+t,\quad \beta_1(t)=a_0+\log D(t),\quad D(t)=1-b(e^t-1),\quad
\alpha_2=a_2,\ \beta_2=b_2.
\]
Thus sources 0, 2, 3 and 4 are identical at every path point, not merely first-order stationary.

The envelope check uses \(|L_\delta|\le A\), \(0<L_\delta'<1\), and exponential saturation of \(L_\delta'\) and \(L_\delta''\) away from \([-A,A]\). Together with
\[
v'(z)^2\le\delta^2+4\delta\sqrt\theta|v(z)|,\qquad |v''(z)|\le2\delta\sqrt\theta,
\]
it yields \(e^{-A}\le b\le e^A\), \(0<u'<1/1000\), \(|u''|<2/3\), and
\[
.8e^x<g_1'(x)<1.25e^x,\quad |g_1''(x)|<2e^x,\quad |g_1'''(x)|<3e^x.
\]
Normalization, positivity, the limit at minus infinity and real analyticity follow from the displayed integral. In particular this is an admissible hard profile, with slack in the original inequalities. Also \(D\in[161/171,19/18]\), \(\alpha_1\in[.7,.8]\), and \(\beta_1\in(11/16,13/16)\); the full segment stays in the original coefficient box.

For source 1, let \(H=h'\), \(d_t(z)=v(z)-u(z+t)\), and \(\lambda_t=be^t/D(t)\in[81/100,200/161]\). Direct differentiation gives
\[
\dot F_t(z)=be^{z+a_0+t}\{e^{u(z+t)}-e^{v(z)}\},\qquad
w_t(z)=\frac{\dot F_t}{F_t'}=\frac{\lambda_t(1-e^{d_t})}{\lambda_t+He^{d_t}},\quad S_t(z)=zw_t(z)-w_t'(z).
\]
The Jacobian derivative \(-w_t'\) is essential. The proof keeps the cancellation in \(|\partial_H w|\le |w|/H\); bounding numerator and denominator terms separately would lose the small central score.

Split at \(r_\delta^2=A/(16\delta)\). On the center, \(|v(z+t)|\le A/2\); the saturation errors are exponentially small. The retained estimates give a central score norm at most \(30\delta(|t|+e^{-A/(2\delta)})\). Globally \(|w_t|\le2\) and \(|S_t(z)|\le3(1+|z|)\); the Gaussian tilted-moment tail estimate gives at most \(9e^{-r_\delta^2/8}\). Consequently, with \(\epsilon_\delta=e^{-1/(1280\delta)}\),
\[
\|S_t\|_{L^2(\phi)}\le40(\epsilon_\delta+\delta|t|).
\]
This controls the whole segment. Its score is strictly positive in norm at \(t=0\) for each positive \(\delta\).

Every response density has common support \((0,\infty)\), with \(p_t(y)=\phi(F_t^{-1}y)/F_t'(F_t^{-1}y)\). To justify integrating the square-root likelihood derivative, the v2 proof changes to the fixed coordinate \(Z_t(z)=F_t^{-1}(F_0(z))\), with Jacobian \(J_t\). The flow equations \(\partial_tZ_t=-w_t(Z_t)\), \(\partial_t\log J_t=-w_t'(Z_t)\) give \(|Z_t-z|\le.1\) and a common integrable envelope proportional to \((1+|z|)^2e^{1+|z|}\phi(z)\). This supplies \(L^2\)-continuous differentiability and DQM, rather than only a formal score. Under \(H^2(P,Q)=\int(\sqrt p-\sqrt q)^2\),
\[
H(P_{1,t},P_{1,0})\le20\epsilon_\delta|t|+10\delta t^2.
\]

The target derivative is
\[
I'(t)=be^{a_0+t}\!\left(Ee^Z[e^{a_2+u(Z+t+a_2)}-e^{u(Z+t)}]-m_h(e^{b_2}-1)\right),\quad m_h=Ee^{h(Z)}.
\]
Oddness of \(h\) gives \(m_h\ge e^{1/2}\). Bounding \(|u|\le A\), the bracket is at most
\(e^{1/2}(e^{.7}-e^{-.1}-e^{.9}+1)\le-.1e^{1/2}\).
Since \(be^{a_0+t}\ge e^{.6}\),
\[
I'(t)\le-k,\qquad k=e^{1.1}/10.
\]
The profile integral is not differentiated as though its internal center moved with \(t\).

For a local perturbation \(t=v/\sqrt N\), source 1 contributes information \((n/N)E S_0^2=E S_0^2/5\le320\epsilon_\delta^2\). Cauchy–Schwarz in the **complete** score space implies
\(V_{\rm eff}\ge k^2/(320\epsilon_\delta^2)\), giving the efficient-variance constant in Section 1. The parametric path is inside the original model, so it is a valid obstruction for the larger score problem.

For the finite lower bound, the other four groups contribute affinity exactly one. Product affinity and \(\mathrm{TV}\le H\) give
\[
\mathrm{TV}(\mathbb P_t,\mathbb P_0)\le\sqrt n(20\epsilon_\delta|t|+10\delta t^2).
\]
Choose
\[
t_N=\min\{1/20,\ (80\sqrt n\epsilon_\delta)^{-1},\ (40\delta\sqrt n)^{-1/2}\}.
\]
Both varying terms are at most \(1/4\), so TV is at most \(1/2\), and target separation is at least \(kt_N\). The nearer-target test has squared error at least one quarter of squared separation when wrong. Maximum risk is at least half the sum, and the sum of test errors is at least \(1-\mathrm{TV}\). Hence maximum risk is at least \(k^2t_N^2/16\). Substituting \(n=N/5\) yields the stronger exact bound
\[
\mathcal R_{N,\delta,\theta}\ge\frac{k^2}{16}\min\left\{\frac1{400},\frac1{1280N\epsilon_\delta^2},\frac{\sqrt5}{40\delta\sqrt N}\right\}.
\]
Factoring the smallest scalar coefficient, \(1/1280\), gives the boxed lower bound. These checks reproduce [lower_bound, Sections 1–6](evidence/docs/interaction_recoverability/warp_family_lower_bound/v2/lower_bound.md), including allocation and loss constants.

## 4. New analytic corollary: polynomial counts separate the profile classes

Let \(N_\delta=5\lceil\delta^{-4}/5\rceil\). For sufficiently small \(\delta\),
\[
\delta^{-4}\le N_\delta<\delta^{-4}+5\le2\delta^{-4},
\quad \frac\delta{\sqrt2}\le\frac1{\delta\sqrt{N_\delta}}\le\delta,
\quad \frac{\delta^2}{2}\le\frac1{N_\delta\delta^2}\le\delta^2.
\]
Here the last inequality in the first line holds for \(\delta\le5^{-1/4}\). To make “sufficiently small” explicit without numerical evaluation, use \(e^x\ge x^5/5!\): if \(\delta\le(240\cdot640^5)^{-1}\), then
\[
e^{1/(640\delta)}\ge\frac{\delta^{-5}}{120\cdot640^5}\ge2\delta^{-4}\ge N_\delta.
\]
For \(0<\delta\le\min\{\delta_*,1,5^{-1/4},(240\cdot640^5)^{-1}\}\), both other terms of the minimum are at least one and its third term lies between \(\delta/\sqrt2\) and \(\delta\). Therefore, uniformly in \(\theta\),
\[
\boxed{\mathcal R_{N_\delta,\delta,\theta}\ge(C_*/\sqrt2)\delta,\qquad
(N_\delta\delta^2)^{-1}=\Theta(\delta^2).}
\]
The deliberately conservative threshold is sufficient, not optimal.

The historical [problem](evidence/docs/interaction_recoverability/problem.md) and [proofs, restricted minimax-rate argument](evidence/docs/interaction_recoverability/proofs.md) concern the special **known** pair \(g_1=g_2=e^x\), \(\theta=1\), and \(0<\delta\le.1\). Their two-quantile inverse has determinant proportional to \(\sinh(\delta s_1(r))\), inverse sensitivity \(O(\delta^{-1})\), and an integrated squared-error argument, giving risk order \(\min\{1,(n\delta^2)^{-1}\}\). With balanced total size this is \(\Theta(\min\{1,(N\delta^2)^{-1}\})\); allocation changes constants only.

If a \(\delta\)-independent finite constant \(C\) gave that upper order over the full original class, at \(\theta=1,N=N_\delta\) it would imply \(\mathcal R\le C\delta^2\). The new inequality contradicts it for \(\delta<C_*/(\sqrt2 C)\) within the displayed range. Thus the historical exponential-profile upper rate cannot extend uniformly over the full original class. This does not contradict the historical subclass theorem or any fixed-\(\delta\) pointwise root-\(N\) statement.

The hard profile pair depends on known \((\delta,\theta)\) but not on \(N\), and is identical at both endpoints and along their path. **It could even be disclosed to the estimator.** Such identical side information changes neither likelihoods nor the testing inequality. The corollary therefore separates geometries within the full profile class; it does **not** isolate an additional cost caused specifically by learning unknown profiles.

## 5. Dependency and defect assessment

| Claim or check | Exact preserved location | Assessment and affected scope |
|---|---|---|
| Full truth class, five laws, score convention | [Family problem](evidence/docs/interaction_recoverability/warp_family/v1/problem.md) | All original coefficient/profile restrictions retained; known warp only; admissible boundary paths qualified. |
| Observable compact inversion and calibration | [Commutator](evidence/docs/interaction_recoverability/warp_family/v1/commutator_and_calibration.md), §§1–6, (1)–(11) | Six-term cancellation, contraction, normalized transports and fixed charts are actual dependencies; central inversion alone is insufficient. |
| Pointwise expansion and arithmetic | [Attainment](evidence/docs/interaction_recoverability/warp_family/v1/attainment.md), §§1–8 | No fatal defect found. Accumulated remainders, moving queries and arithmetic assumptions deserve external proof scrutiny. No empirical backend or uniform risk follows. |
| Complete variance upper | [Conditioning](evidence/docs/interaction_recoverability/warp_family/v1/conditioning.md), §§2–5 | Includes remote calibration and coefficient influence; using only central variance would invalidate the claimed full bound. |
| Hard profile and nonlinear path | [Construction](evidence/docs/interaction_recoverability/warp_family_lower_bound/v2/construction.md), (1)–(6) | Profiles really remain fixed; moving their internal center would be a different, invalid argument. |
| Global envelopes, actual score, DQM, margin | [Admissibility](evidence/docs/interaction_recoverability/warp_family_lower_bound/v2/admissibility_and_score.md), §§1–7 | Checked global tails and common-coordinate domination, not just a formal derivative at zero. No fatal defect found. |
| Five-source information and finite risk | [Lower bound](evidence/docs/interaction_recoverability/warp_family_lower_bound/v2/lower_bound.md), §§1–6 | Allocation \(1/5\), product affinity and testing factor \(1/16\) retained. Dropping the quadratic path term would overclaim. |
| Corollary | This document, §4 | New analytic deduction from the finite inequality; independent of a triangular-array CLT. |
| Earlier review packet | [v2 brief](evidence/docs/interaction_recoverability/external_review/v2/review_brief.md), [errata](evidence/docs/interaction_recoverability/external_review/v2/errata.md) | Historical compact/local review material, not external validation of this family theorem. Kept byte for byte. |
| Reference efficiency and feasibility | [Efficiency](evidence/docs/interaction_recoverability/efficiency_audit/v1/efficiency_argument.md), [no-go](evidence/docs/interaction_recoverability/practical_feasibility/v1/feasibility_decision.md) | The old \(\delta=.1\) efficient variance and efficient attainment remain OPEN; reported variance diagnostics are not certified errors. Backend remains NO-GO. |

Concrete objections to a broader contribution claim are the extremely small sufficient separation range; known Gaussian latent law and known warp; strong calibrated, labeled private anchors; no finite-sample upper risk, useful count regime or efficient estimator; no biological calibration; and no new theorem-covered learned GNN. These limit significance without refuting the stated mathematical example. The score/tail construction and global inversion need independent checking before publication.

The crossover \(N\asymp\delta^2e^{1/(320\delta)}\) is obtained by equating two conservative lower-bound terms. Neither that location nor the individual regimes are proved minimax rates or a CLT onset. At fixed desired accuracy, once the quadratic-path term falls below a threshold this lower inequality ceases to exclude that accuracy, even if its information term is large. An exponential variance constant therefore does not prove exponential sample requirements for fixed accuracy.

## 6. Exact prior overlap, with versions inspected

Primary text was read from the retained local caches first. Bounded HTML retrieval supplemented Kaji, Heinrich–Kahn, DExtrI and Yao–de la Llave. No required full text was unavailable. The packet does not redistribute those papers; it supplies stable version links. Statements about relevance below are this review's deductions from the hypotheses, not claims made by the cited authors.

| Primary source and inspected result | Observation model, assumptions and conclusion | Overlap and limit of transfer |
|---|---|---|
| [Trabs, arXiv:1307.6610v2](https://arxiv.org/abs/1307.6610v2), Theorem 2.7 | Regular indirect LAN experiments; differentiable target; score operator and its adjoint. Adjoint range characterizes regular estimability, with a minimum-norm information bound and convolution result. | Supplies the efficiency language and information inequality. It does not construct this empirical statistic, control its growing query scheme, or give the present five-law finite-risk bound. |
| [Kaji, arXiv:1910.07572v1](https://arxiv.org/html/1910.07572v1), Theorem 4.1, Proposition 4.2, Theorem 5.1 | Inverse-map differentiability for distributions smooth away from finitely many jumps with positive density; transformed quantiles require smooth transforms and tail integrability. The integrable empirical/quantile-process result and fixed-class L-statistic calculus permit specified random integrators. | Provides established quantile delta-method machinery. The five-source calibration, changing queries, growing series and explicit tail/remainder schedules still need the retained separate proof. It is not a ready-made uniform \(\delta\downarrow0\) MSE theorem. |
| [Heinrich and Kahn, arXiv:1507.04313v1](https://arxiv.org/html/1507.04313v1), Theorems 3.2, 3.3, 3.5 | Finite mixtures near a distribution with fewer components. Smoothness/nonzero-derivative assumptions yield local lower bounds; strong identifiability and stronger smoothness yield matching local/global estimators. Risk is Wasserstein loss for the mixing distribution. Pointwise root-\(n\) estimation is nonuniform near collisions. | Establishes the general distinction between pointwise regularity and uniform difficulty near singularities. The mixture parameter, loss, class and exponents differ; no result there implies this exponential variance constant or the present scalar-risk bound. Version v1 was inspected; no claim is made about changes in later versions. |
| [DExtrI, arXiv:2608.19849v1](https://arxiv.org/html/2608.19849v1), Theorems 3.3–3.4 and Appendix D | Shared-noise conditional laws on nondegenerate coordinate-axis intervals. Main assumptions include analytic/nonpolynomial profiles, normalized monotone warps, support/nonproportionality, pooled growth and real-frequency exponential-polynomials. Conclusions identify parameterizations modulo indeterminacies and determine off-axis population laws. | Close structural antecedent; five isolated calibrated laws and scalar statistical recovery differ from continuum identification. Appendix D.3 permits linear independence plus pair rigidity; D.8 gives a sparse-foliation rigidity route; D.9 uses complex singularities of profile derivatives. Exclusion from the main theorem alone is insufficient for novelty. Our envelope does not imply those rigidity conditions; warp branch points do not establish profile-derivative singularities, and exponential derivatives are entire. D.12 treats bounded-degree polynomial profiles under pooled finite-power independence; our envelope excludes polynomials. D.13 treats affine warps under F1–F5, identifying directions/supports and profile classes modulo polynomials; our second warp is nonaffine. None of these inspected results supplies this finite-five-source estimator or finite-risk inequality. Exact priority remains open. |
| [Bennett et al., arXiv:2208.08291v3](https://arxiv.org/abs/2208.08291v3), Theorem 2 and adjacent identification discussion | Linear functionals of solutions to conditional moment equations. Identification uses a closure condition; debiasing uses an adjoint-range condition, while strong identification imposes the stronger \(\mathrm{Range}(P^*P)\) condition. | Full-profile inverse stability is not required for every target functional. This established principle is not a new contribution here. The paper's operator and data model are not the current five-law density experiment. |
| [Yao and de la Llave, arXiv:2110.15893v1](https://arxiv.org/html/2110.15893v1), equations (15)–(17), §3.4 | Weighted cohomological equations solved by iterates/Neumann series under contraction conditions, with accelerated summation. | Iterated inversion is established machinery. Here the normalized difference quotient uses composition by \(C\), contracting a centered Hölder seminorm by \(\lambda^\gamma\); the multiplier is one, so a theorem requiring multiplier sup-norm below one cannot be quoted verbatim. The calibrated observable identity and its empirical remainders are the specific work. |

This comparison supports a concrete calibrated example, not a new general adjoint-range, quantile, cohomological, Riesz-learning or testing principle. DExtrI's broader appendix is materially relevant and prevents an argument based only on its main theorem's narrower assumptions. Establishing how much significance the finite-design construction adds beyond these tools is an external-review question.

## 7. Open external questions

The [unsent request](reviewer_request.md) asks an expert to verify the full upper proof and hard path, check exact overlap including DExtrI's appendix, and assess the value of the calibrated example. Sharp risk rates, leading exponential constants, a useful uniform shrinking-separation upper theorem, efficient estimation, the old reference-efficiency problem, finite-sample inference, fixed-resolution realization and biological applicability remain open. This review supplies no new research or fitting cycle.
