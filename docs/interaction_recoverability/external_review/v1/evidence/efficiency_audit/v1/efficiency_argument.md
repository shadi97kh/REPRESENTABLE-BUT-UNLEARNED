> Package note — copied historical evidence. Mathematical prose and formulas are retained; links to excluded repository artifacts are rendered as inactive references, and any link changes are recorded in the package manifest. This is not a new proof or an edit to the original file.

# Bounded efficiency audit of the complete five-law experiment

**Decision: unresolved efficiency; the existing backend no-go remains.** This audit proves a full-model-compatible known-median correction and genuine nested-submodel lower bounds. The deterministic numbers leave a gap of about 38,776-fold. Removing this one known moment does not establish the efficient influence, an efficient estimator, or a substantive new learning contribution.

## Experiment, score and target

**PROVED HERE / saved model retained.** Use the [original model](../../unknown_profiles/model.md), [analytic tangent specification](../../local_five_law/v1/proof.md), and [compact construction](../../compact_attainment/v1/construction.md) without restriction to the numerical dictionary. There are five independent response groups of size `n`, `N=5n`, with weights \(\pi_a=1/5\). The known warp is \(h(z)=z+.1z\sqrt{1+z^2}\), \(Z\sim N(0,1)\). Both normalized analytic profiles and all four coefficients remain unknown. Only for reference geometry,
\[
g_1=g_2=e^x,\qquad
\alpha=(2/3,3/4),\quad \beta=(3/5,4/5).
\]
Write the five shift pairs as
\[
(s_a,t_a)=(0,0),(\alpha_1,\beta_1),(\alpha_2,\beta_2),(1,0),(0,1),
\quad F_a(z)=g_1(z+s_a)+g_2(hz+t_a).
\]
A full tangent is \(v=(r_1,r_2,d_{\alpha_1},d_{\beta_1},d_{\alpha_2},d_{\beta_2})\).
Put \(d_{s,a}=d_{\alpha_a},d_{t,a}=d_{\beta_a}\) for \(a=1,2\), and zero otherwise. Its response derivative is
\[
q_a(z)=r_1(z+s_a)+r_2(hz+t_a)
 +g'_1(z+s_a)d_{s,a}+g'_2(hz+t_a)d_{t,a}.
\]
At the reference,
\[
\begin{split}
q'_a&=r'_1(z+s_a)+h'r'_2(hz+t_a)
 +e^{z+s_a}d_{s,a}+h'e^{hz+t_a}d_{t,a},\\
F'_a&=e^{z+s_a}+h'e^{hz+t_a},\\
F''_a&=e^{z+s_a}+\{h''+(h')^2\}e^{hz+t_a},\\
h''(z)&=.1z(3+2z^2)/(1+z^2)^{3/2}.
\end{split}
\]
Differentiating the density at a fixed observed outcome, rather than at a moving latent quantile, gives
\[
\boxed{(Sv)_a(F_a(z))
=-\frac{(\phi q_a/F'_a)'(z)}{\phi(z)}
=z\frac{q_a}{F'_a}-\frac{q'_a}{F'_a}
+\frac{q_aF''_a}{(F'_a)^2}.}\tag{1}
\]
These centered scores have finite Gaussian \(L^2\) norms for the original bounded-weight analytic directions. The saved fixed-path DQM argument applies to the two-sided paths below. The experiment's tangent space is
\(\mathcal T=\overline{\operatorname{span}\{Sv:v\ {\rm admissible}\}}\)
in \(\mathcal H=\bigoplus_{a=0}^4 L^2_0(P_a)\), with
\(\langle f,g\rangle_\pi=\frac15\sum_aE_af_ag_a\).
There is no finite-sieve substitution.

Let \(\Delta_{u,v}r(x)=r(x+u+v)-r(x+u)-r(x+v)+r(x)\). The full target derivative is
\[
DI[v]=E\Delta_{\alpha_1,\alpha_2}r_1(Z)
 +E\Delta_{\beta_1,\beta_2}r_2(hZ)+b_c^\top d_c,\tag{2}
\]
where the coefficient order throughout this audit is
\((\alpha_1,\beta_1,\alpha_2,\beta_2)\). With
\(m_0=e^{1/2}\), \(m_h=Ee^{hZ}\),
\[
b_c=\bigl(m_0e^{\alpha_1}(e^{\alpha_2}-1),
m_he^{\beta_1}(e^{\beta_2}-1),
m_0e^{\alpha_2}(e^{\alpha_1}-1),
m_he^{\beta_2}(e^{\beta_1}-1)\bigr).\tag{3}
\]
For an independent formula for the profile entries, put
\(M_1(w)=e^{w^2/2}\), \(M_2(w)=Ee^{whZ}\), and
\(D_j(w)=M_j(w)(e^{wc_{j1}}-1)(e^{wc_{j2}}-1)\).
The sine entry is \(\Im D_j(1+ik)\); the cosine-minus-one entry is
\(\Re D_j(1+ik)-D_j(1)\). All integrals exist, since the positive-tail exponent is at most \(-.4z^2+O(z)\).

## Full-compatible reference and efficiency

**STANDARD RESULT APPLIED.** Define
\[
V_{\rm eff}=\inf\{\|\psi\|_\pi^2:
 \langle\psi,Sv\rangle_\pi=DI[v]\text{ for every admissible }v\}.\tag{4}
\]
The definition concerns regular local root-\(N\) inference. The closed affine solution set is nonempty by the saved construction; its minimum is \(P_{\mathcal T}\psi\). It does not assert existence of a feasible efficient estimator or a finite-sample MSE bound. The adjoint-range/minimum-norm formulation is established theory: see [Trabs, Theorem 2.7](https://arxiv.org/html/1307.6610v2), reread for this audit.

**PROVED HERE — compatibility audit of the saved operator.** The compact representer uses
\[
\kappa_{a,z}(y)=5\,\frac{F'_a(z)}{\phi(z)}
 \{\Phi(z)-\mathbf1[y\le F_a(z)]\}
\]
in coordinate \(a\), zero in the others. Integration of (1) gives
\(\langle\kappa_{a,z},Sv\rangle_\pi=q_a(z)\).
The factor 5 is essential. Normalization gives
\(r_1(-1)=q_0(-1)-q_3(-1)\). The compact 18-term telescope satisfies
\(E_{45}(q;x)=r_1(Cx)-r_1(x)\), \(C=h^{-1}(h+5)-4\).
Its convergent difference series, with the saved transports, supplies profile evaluation duals \(u_x,v_y\), where \(v_y=\kappa_{0,h^{-1}y}-u_{h^{-1}y}\).

The actual residuals and complete influence are
\[
\zeta_{i,\pm}=\kappa_{i,\pm1}
 -\kappa_{0,h^{-1}(\pm h(1)+\beta_i)}
 +u_{h^{-1}(\pm h(1)+\beta_i)}-u_{\pm1+\alpha_i},
\]
\[
\chi_i=J_i^{-1}(\zeta_{i,-},\zeta_{i,+})^\top,\qquad
J_i=\begin{pmatrix}
e^{-1+\alpha_i}&e^{-h(1)+\beta_i}\\
e^{1+\alpha_i}&e^{h(1)+\beta_i}
\end{pmatrix},\qquad
\psi=\psi_A+b_c^\top\chi.\tag{5}
\]
Here \(\psi_A\) is the compact outward target representer in the linked construction and [theorem](../../compact_attainment/v1/theorem.md).
The common profile normalizer and central reference cancel in each residual difference. The anchor coordinates have no coefficient score; the residuals pair to \(J_i d_i\) and vanish against every pure profile direction. Thus (5) satisfies all profile and coefficient equations in (4). Each source coordinate has finite individual \(L^2\) norm by the saved compact-series and outward-tail bounds. Continuity extends compatibility to the full closure.

The inspected covariance code (omitted repository reference: `../../practical_feasibility/v1/influence_variance.py`) forms all signed contributions within a source before squaring, transforms both residual pairs by \(J_i^{-1}\), and then applies (3). It includes source-zero covariance between the anchor and all four coefficients. Its covariance factor is \((1/5)5^2=5\). No sign, ordering, or omitted shared-source term was found in these inspected operations. This is an operator audit, not a new verification of every earlier stochastic attainment remainder.

Consequently \(V_{\rm eff}\le V_{\rm constructed}:=\|\psi\|_\pi^2\). Interpreting the previously reported large variance as an unavoidable statistical lower bound would reverse this inequality.

## Exact known-median projection

**PROVED HERE.** Every admissible model has
\(F_0(0)=g_1(0)+g_2(0)=2\). Monotonicity gives
\(P_0(Y\le2)=P(Z\le0)=1/2\). Define
\[
Q_0(y)=\mathbf1[y\le2]-1/2,\qquad Q_a=0\ (a\ne0).
\]
It is centered and \(\|Q\|_\pi^2=(1/5)(1/4)=1/20\).
Along every admissible DQM path, \(E_{0,t}Q_0=0\) identically. Differentiating this bounded expectation gives
\(\langle Q,Sv\rangle_\pi=0\), including all four coefficient directions. The bounded linear pairing extends to \(\mathcal T\); hence \(Q\in\mathcal T^\perp\).

For \(c=\langle\psi,Q\rangle_\pi\), orthogonal projection gives
\[
\lambda_*=20c,\qquad \psi_c=\psi-\lambda_*Q,\qquad
U:=\|\psi_c\|_\pi^2=V_{\rm constructed}-20c^2.\tag{6}
\]
It follows exactly that \(V_{\rm eff}\le U\le V_{\rm constructed}\).
Strict reduction requires \(c\ne0\). The computed covariance is decisively nonzero as a diagnostic; complete numerical certification of that fact is not asserted.

For a raw source-zero signed atom \(w\kappa_{0,z}\), its contribution to \(c\) is
\[
-\frac{w}{2}\frac{F'_0(z)}{\phi(z)}\Phi(-|z|).\tag{7}
\]
Indeed \(E[(p-\mathbf1[U\le p])(\mathbf1[U\le1/2]-1/2)]
=-(\min(p,1/2)-p/2)\).
The source weight cancels the kernel's factor 5. Formula (7) is accumulated for the complete anchor/residual vector before the same coefficient transformation and target gradient as (5). It does not compute only the anchor covariance.

**PROVED HERE — observable adjustment, distinct from an efficient estimator.** For a fixed public scalar \(\lambda\), the adjustment to the saved mathematical statistic is
\[
-\frac{\lambda}{N}\sum_{i=1}^nQ_0(Y_{0i})
=-\frac{\lambda}{5}\overline Q_0.\tag{8}
\]
Using \(-\lambda\overline Q_0\) would be wrong by a factor of five. This uses only observed source-zero responses and the known normalization. Its influence is \(\psi_\eta-\lambda Q\), with variance
\(V_\eta-2\lambda c_\eta+\lambda^2/20\).
A reference number declared as a public fixed constant is implementable; it need not improve variance at other truths. The optimal reference \(\lambda_*\) in (6) is a geometry quantity, not an unknown truth supplied to a learner. At the reference, a frozen rounded value \(\lambda_0\) has exact variance \(U+(\lambda_0-20c)^2/20\), not exactly \(U\) unless the two weights agree.

If a source-based estimate satisfies \(\widehat\lambda\to_P\lambda_\eta\), its extra remainder is
\[
(\widehat\lambda-\lambda_\eta)N^{-1}\sum_iQ_0(Y_{0i})
=o_P(N^{-1/2}).
\]
Independence or cross-fitting is not needed for this product once consistency is proved. **OPEN:** the saved artifacts do not prove source-only estimation of the optimal full-influence covariance weight throughout the unknown-profile class. No such fit, backend or efficient estimator was implemented here. The exact median moment survives all admissible local alternatives; a pointwise consistency assertion alone should not silently be upgraded to uniform inference.

## Frozen submodels and admissibility

**PROVED HERE.** The dictionary was frozen in calculate.py (omitted repository reference: `calculate.py`) and its ledger hash before any new geometry output. Order:
\[
\alpha_1,\beta_1,\alpha_2,\beta_2;\quad
(g_1:\sin x,\cos x-1;\ g_2:\sin x,\cos x-1)e^x;
\quad(g_1:\sin2x,\cos2x-1;\ g_2:\sin2x,\cos2x-1)e^x.
\]
The nested sizes are exactly 4, 8 and 12. No direction was added or selected using results.

For derivative order \(j\le3\), \(k=1,2\),
\[
|\partial^j(e^x\sin kx)|\le(1+k^2)^{j/2}e^x,\quad
|\partial^j(e^x(\cos kx-1))|
\le\{(1+k^2)^{j/2}+1\}e^x<13e^x.
\]
Every direction is real analytic, vanishes at zero, and tends to zero at negative infinity. A joint parameter cube with all amplitudes strictly smaller than \(1/5200\) keeps each profile's four-direction weighted-\(C^3\) deviation below \(4(13)/5200=.01\). In particular it preserves positivity, the original derivative envelopes, and normalization. Coefficient perturbations of this size remain strictly inside the declared coefficient neighborhood. This proves a genuine two-sided finite-dimensional submodel through the reference, rather than merely formal tangent expressions.

For each dictionary let
\[
G_{jk}=\frac15\sum_{a=0}^4\int(Sv_j)_a(F_a(z))(Sv_k)_a(F_a(z))\phi(z)\,dz,\qquad
b_j=DI[v_j].\tag{9}
\]
All five laws enter each profile-direction Gram; only sources 1 and 2 respond to coefficient changes. The code implements (1) with both derivative terms and these exact supports.

**PROVED HERE — exact rank, independently of numerical tolerance.** In fact each of these Grams is positive definite. If a dictionary combination has zero score, (1) implies
\((\phi q_a/F'_a)'=0\). The profile and coefficient envelopes make \(\phi q_a/F'_a\to0\) in each tail, so \(q_a=0\) for every source. The anchor telescope then gives \(r_1(Cx)=r_1(x)\) on the compact invariant interval. Contraction and continuity make \(r_1\) constant there. Analyticity extends constancy to the real line; \(r_1(0)=0\) makes it zero. Source zero forces \(r_2=0\). For sources 1 and 2 the remaining equation is
\(e^{z+\alpha_i}d_{\alpha_i}+e^{hz+\beta_i}d_{\beta_i}=0\).
Nonconstant \(h(z)-z\) forces both coefficients to vanish. Finally the four trigonometric functions within each profile are linearly independent. Thus the original dictionary coefficients all vanish.

This proves no exact zero-information direction occurs in these finite submodels. It does not establish a lower spectral bound for an infinite profile dictionary.

**STANDARD RESULT APPLIED.** For a general PSD Gram, full influence compatibility implies \(b\perp\ker G\): a null combination has zero score and therefore zero pairing with \(\psi\). The minimum-norm representer on its finite score span is
\(r_d=\sum_j(G^\dagger b)_jSv_j\), with squared norm
\[
L_d=b^\top G^\dagger b.
\]
For the present Grams \(G^\dagger=G^{-1}\). Every full-compatible influence projects to \(r_d\), so
\[
\boxed{L_4\le L_8\le L_{12}\le V_{\rm eff}\le U\le V_{\rm constructed}.}\tag{10}
\]
These are exact inequalities between the defined population quantities. A twelve-equation solution is only a restricted-submodel representer and supplies no upper bound for the full model. No ridge-regularized inverse is substituted for (9).

## Deterministic values and numerical qualifications

**NUMERICALLY CHECKED — approximations, not certified endpoints.** raw.json (omitted repository reference: `raw.json`) retains both frozen configurations, all five source Gram matrices, complete target derivatives, eigenvalues, coefficients, residuals and covariance components. The original full variance is read from the preserved final variance record (omitted repository reference: `../../practical_feasibility/v1/variance_final.json`); it is not rerun.

| Quantity | Approximation |
|---|---:|
| Constructed variance \(V\) | 79,390,122.7930676 |
| Anchor-only median covariance | -0.0862098813 |
| Four-coefficient correction's median covariance | 56.8894759144 |
| Complete \(c\) | 56.8032660332 |
| Reference \(\lambda_*=20c\) | 1,136.0653206631 |
| Removed variance \(20c^2\) | 64,532.2206407 |
| Removed fraction | 0.08128495% |
| Corrected variance \(U\) | 79,325,590.5724269 |
| Corrected standard-deviation constant | 8,906.4914850 |
| \(L_4\) | 1,568.5179383034 |
| \(L_8\) | 1,981.9138167828 |
| \(L_{12}\) | 2,045.7389998777 |
| \(U-L_{12}\) | 79,323,544.8334 |
| \(U/L_{12}\) | 38,776.0074 |

Using only the anchor would give the wrong sign and miss almost all the median covariance. Both new configurations give the same practical conclusion. The new covariance quadrature splits additionally at zero; the reused variance quadrature did not. They approximate the same continuously integrated truncation, but the reported subtraction is not literally the squared norm of one common finite atom collection.

The score integrals use composite Gauss rules on \([-12,12]\): panel width/order .25/8 and .125/12. The median operator uses fixed series depth 512 and lattice radius 24, central orders 128/256 and outer panel widths .125/.0625. Outer integration splits at the known median and the outward-weight join. There is no expansion beyond these two predetermined configurations.

| Dictionary size | Smallest numerical eigenvalue | Condition number | Numerical rank | Absolute rank threshold |
|---|---:|---:|---:|---:|
| 4 | .001402827645 | 60.3698 | 4 | \(8.469\,10^{-14}\) |
| 8 | .001154599321 | 487.9263 | 8 | \(5.634\,10^{-13}\) |
| 12 | .001017645762 | 1599.3707 | 12 | \(1.628\,10^{-12}\) |

The threshold is \(10^{-12}\lambda_{\max}\). No eigenmode was dropped; the reported primary value uses an unregularized solve. The thresholded pseudoinverse diagnostic agrees within approximately \(5\,10^{-11}\). Maximum solve residual is \(7.24\,10^{-15}\); discarded-\(b\) norm is zero. These numerical ranks agree with the analytic positive-definiteness argument; they do not prove it.

Refinement changes the Gram in operator norm by \(4.628\,10^{-16}\), \(b\) by \(1.231\,10^{-15}\), and the three bounds by at most \(3.616\,10^{-11}\). Complete median covariance changes by \(2.092\,10^{-8}\). The largest integrated score mean is \(1.637\,10^{-16}\), and the direct versus complex-moment target check differs by \(1.333\,10^{-15}\). The newly integrated coefficient gradient differs from the preserved covariance calculation's gradient by \(1.937\,10^{-13}\). All are diagnostics, not rigorous integration or floating-point errors. Extended accumulation still uses float64 CDFs and initial Gauss nodes.

**PROVED HERE — separate Gram tail control.** For each individual frozen direction, \(|q/F'|\le2\), \(|q'/F'|\le4\), and
\(F''/F'\le1.3+.2|z|\). Consequently \(|Sv|\le7+3|z|\).
Every omitted Gram entry at radius \(R\) is bounded by
\[
E[(7+3|Z|)^2\mathbf1(|Z|>R)]
=116\Phi(-R)+(84+18R)\phi(R).
\]
At \(R=12\) this is below \(10^{-28}\), using Mills' inequality and elementary exponential bounds. The resulting \(d\)-dimensional Gram operator tail is at most \(d\,10^{-28}\). This controls only the omitted tails; the retained quadrature is not certified. For each target derivative, the absolute four-shift integrand is at most \(32e^z\) or \(32e^{h(z)}\); coefficient factors are smaller. At \(z\ge12\), \(h(z)\le .1z^2+z+.05\), so
\[
\int_{12}^{\infty}e^{h(z)}\phi(z)\,dz
\le e^{-45.55}/8.6<e^{-45}/8.
\]
The negative tail and the first-profile moment tail are smaller. Each omitted target-derivative entry is therefore below \(8e^{-45}<10^{-18}\), and its vector tail norm is below \(\sqrt d\,10^{-18}\). These tail bounds do not certify retained quadrature or roundoff. In particular the profile-2 direct/complex check shares integration nodes: it independently checks algebra, not numerical integration.

**PROVED HERE — propagate the saved influence error correctly.** For the exact continuously integrated truncation \(\psi_T\) at \(K=512,R=12,M=24\), the saved [analytic error analysis](../../practical_feasibility/v1/variance_analysis.md) gives
\(\|\psi-\psi_T\|_\pi\le\epsilon<5.905\).
The central, response and lattice contributions are separately bounded there. Since \(P_Qf=20\langle f,Q\rangle_\pi Q\),
\[
|c-c_T|\le\epsilon/\sqrt{20}<1.320399,\qquad
\|(I-P_Q)\psi-(I-P_Q)\psi_T\|_\pi\le\epsilon.\tag{11}
\]
Thus the *corrected* exact truncation norm has the same 5.905 error bound. Neither \(8,906.4915\pm5.905\) nor the displayed numeric chain is a certified interval: the numerical integration error is still missing.

More explicitly, if numerical certificates established
\(|V_T-\widetilde V|\le\delta_V\) and
\(|c_T-\widetilde c|\le\delta_c\), an upper bound would be
\[
U\le\left[
 \sqrt{\max\{0,\widetilde V+\delta_V
       -20(\max\{0,|\widetilde c|-\delta_c\})^2\}}+\epsilon
 \right]^2.\tag{12}
\]
This propagates both covariance and norm error, including the shared reference source. The current computation does not supply \(\delta_V,\delta_c\). A refinement difference cannot fill that role. In particular strict nonorthogonality would follow from
\(|\widetilde c|>\epsilon/\sqrt{20}+\delta_c\); that numerical certificate was not constructed.

For lower bounds, any fixed rational vector \(t\) obeys
\(L_d\ge2b^\top t-t^\top Gt\).
If certificates supplied \(\|b-\widetilde b\|\le\delta_b\) and
\(\|G-\widetilde G\|_{\rm op}\le\delta_G\), then
\[
L_d\ge2\widetilde b^\top t-t^\top\widetilde Gt
       -2\delta_b\|t\|-\delta_G\|t\|^2.\tag{13}
\]
The recorded solve coefficients can be interpreted as fixed decimal rationals for this purpose. No numerical-rank assumption is needed in (13). The present retained-quadrature errors are not certified, so the table reports approximations to the exact lower bounds, not literal certified numerical lower endpoints.

## Exact unresolved question and contribution

**PROVED HERE / STANDARD RESULT APPLIED.** Since \(Q\perp\mathcal T\),
\(P_{\mathcal T}\psi_c=P_{\mathcal T}\psi=\psi_{\rm eff}\).
The unresolved gap decomposes exactly as
\[
U-L_{12}
=\underbrace{\|\psi_c-P_{\mathcal T}\psi_c\|_\pi^2}_{U-V_{\rm eff}}
+\underbrace{\|P_{\mathcal T}\psi_c-r_{12}\|_\pi^2}_{V_{\rm eff}-L_{12}}.\tag{14}
\]
Neither term has been controlled sharply. This is the specific missing full-closure projection, not a missing finite Gram inverse. A large gap does not say which endpoint is close to efficiency.

The median correction is ordinary orthogonal projection onto a known zero-mean moment. The submodel calculation is ordinary information/Riesz projection onto a finite score span. These established operations are useful for diagnosing this model-specific construction; they are not a new learning algorithm. The saved [pilot comparison](../../learned_target_v1/novelty.md) already identifies regularized finite-score Riesz/one-step correction and residual Galerkin enrichment. Its negative decision (omitted repository reference: `../../learned_target_v1/decision.md`) reports 3.5569% MSE improvement, below the frozen 10% screen. This audit neither reruns that pilot nor repairs its full-class uncertainty and nonlinear-remainder gaps.

**OPEN.** Efficient variance, an attaining efficient estimator, useful finite-sample accuracy, MSE, certified numerical endpoints and biological siRNA validity remain unestablished. The new values leave the [practical backend no-go](../../practical_feasibility/v1/feasibility_decision.md) unchanged. The corrected asymptotic scale remains about 8.906 at \(N=10^6\), without demonstrated CLT onset; this is not a finite-sample standard error or a necessary sample size. No new training or backend cycle follows.

The mathematical/code checks were performed by the primary agent and two separate agents within the same system. This is an internal adversarial cross-check, not an external independent review. Frozen results were preserved; see commands and resources (omitted repository reference: `commands_and_resources.md`).

