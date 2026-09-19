# Static arithmetic audit and local integration certificates

**PROVED HERE / qualifications stated below.** Static inspection supports the saved compact construction's \(N\,\mathrm{polylog}N\) arithmetic-existence bound under its specified input and primitive model. Its diagnostic implementation is not a complete certified backend. In particular, exact classification of arbitrary real inputs at the response-cap boundary is not supplied by a finite rounded encoding. The new [numerical stability argument](numerical_stability.md) treats classification in probability; it does not retroactively make exact boundary classification a finite-bit operation.

This is a static self-review with a separate agent cross-check inside the same assistant workflow, not an external independent review. The saved theorem, construction, algorithm, quadrature proof, code, checks, public count records, prior static review and accounting were read. No numerical computation, sampled input, fit or benchmark was executed for this addendum.

## 1. Dependency audit

**PROVED HERE, by checking the displayed recurrences:** the saved derivative and knot dependencies have the following orders. Constants below can be large but are independent of \(N\).

| Object | Actual dependence and qualification |
|---|---|
| Warp and inverse derivatives through the needed order | The displayed \(h^{(2)},\ldots,h^{(5)}\) bounds are constant; \(h'\) grows linearly in the argument bound. The inverse recurrence uses \(h'\ge1\), so there is no unknown denominator. |
| Derivatives of \(C^k\) on \(J\) | The saved \(D_r\) recursion is uniform in \(k\); it follows from \(u_{r,k+1}\le\lambda u_{r,k}+\sum_{m\ge2}b_m{\cal B}_{r,m}(u_{1,k},\ldots)\). The stronger decaying bounds below are also available. |
| Known weight derivatives | Their public bounds are polynomials in \(R+M+4\), multiplied by finite sums of length \(O(M)\). With \(R,M=O(\sqrt{\log N})\), these remain polynomials in \(\log N\). The unknown profiles are not differentiated in this calculation. |
| Smooth-piece fourth derivatives | The factor \(N^3\) is the worst-case empirical polygon slope \(n\max_i(Y_{(i+1)}-Y_{(i)})\le nN^2\). The remaining displayed factor is polynomial in \(K,R,M\). It is not a measured slope or runtime. |
| Knot count | Each central probability map has range length at most \(2\lambda^k\), hence at most \(2n\lambda^k+2\) crossed knots. Summing the 18 terms yields \(O(N+K)\), not a quadratic rank-pair construction. |
| Knot conditioning | The lower derivative \(m_K=2^{-K-11}\) is small but only polynomially small in \(N\), since \(K=O(\log N)\). Forward bisection therefore requires \(O(\log N)\) precision and iterations at the saved polynomially small gap width. |
| Touching/overlapping knot enclosures | Taking their union and excluding it is valid. The total excluded length is at most the sum of enclosure lengths, even when endpoints coincide. Exact transcendental equality or exact ordering of equal roots is unnecessary. |
| Primitive and summation precision | Known arguments have magnitude polynomial in \(\log N\). Straight-line sensitivity and the number of additions introduce polynomial factors in \(N\); their logarithms are \(O(\log N)\). A final interval enclosure must meet each allocated absolute-error share. The saved code does not implement that enclosure graph. |
| Coefficient iteration | The deterministic amplification \(A_T\) may have \(\log A_T=O((\log N)^2)\). Retaining this factor makes the old sufficient fractional-input/working-precision order consistent. Removing it requires the separate stochastic argument, not an assertion of empirical contraction. |
| Reading and sorting input | On the capped input model, integer bits are \(O(\log N)\). Reading \(O((\log N)^2)\) fractional bits for each of \(N\) heights remains \(N\,\mathrm{polylog}N\) work. An external fixed-resolution recording does not supply these digits. |
| Schedule boundaries | Rational-power tests and outward enclosures avoid exact logarithm-ceiling tests. The already controlled extra lattice shell covers outward integer selection. |

**Concrete qualification — exact cap classification.** The mathematical gate tests the original response against \(0\) and \(N^2\). At an arbitrary real number arbitrarily close to either boundary, a prescribed finite encoding need not determine that comparison. Thus the saved deterministic comparison is conditional on original cap membership/known classification, or on a promise of valid capped input. Its note about a guard band identifies a real additional obligation. The high-probability rounded gate analyzed in this continuation is sufficient for an in-probability estimator comparison, while exact agreement on every input remains a stronger assertion.

**OPEN / implementation boundary.** The delivered `mpmath` code uses approximate comparisons, ordinary floating-point ceilings, ordinary summation, and diagnostic bisection. It contains no complete knot compiler, certified cap classifier, interval primitive library, or end-to-end integration/error allocator. Its unconditional sampled-data refusal is appropriate. No fatal \(N\)-dependence error was found in the conditional arithmetic-existence calculation; no complete executable certificate has been verified.

## 2. Derivatives decay along the contraction

**PROVED HERE.** Let \(\lambda=271/290\). Write \(b_m\) for the saved constant upper bounds on \(|C^{(m)}|\), \(m=2,3,4\); use the saved bounds for the derivatives of the outer maps \(T_j\) as well. Define
\[
E_1=1,\qquad
E_r={1\over\lambda(1-\lambda)}
\sum_{m=2}^{r}b_m{\cal B}_{r,m}
(E_1,\ldots,E_{r-m+1}),\quad r=2,3,4.
\]
Then
\[
\sup_{x\in J}|(C^k)^{(r)}(x)|\le E_r\lambda^k,\qquad
k\ge0,\quad r=1,\ldots,4. \tag{A1}
\]
For \(r=1\) this is contraction. Suppose the lower derivative orders have this form. Faà di Bruno and the nonnegative coefficients of the Bell polynomials give
\[
u_{r,k+1}\le\lambda u_{r,k}+S_r\lambda^{2k},
\quad
S_r=\sum_{m=2}^{r}b_m{\cal B}_{r,m}(E_1,\ldots).
\]
Since \(u_{r,0}=0\) for \(r\ge2\), summing the geometric recursion yields
\(u_{r,k}\le S_r\lambda^{k-1}/(1-\lambda)=E_r\lambda^k\).
This proves (A1) by induction in \(r\). No fitted derivative or empirical contraction is used.

For each central query \(f_{\ell k}(t)=T_j(C^kt)-r\), or the identity outer query when present, compose (A1) with the known outer map and \(\Phi\). There are fixed constants \(A_r\) such that
\[
\sup_{t\in J}|(\Phi\circ f_{\ell k})^{(r)}(t)|
\le A_r\lambda^k,\quad r=1,\ldots,4. \tag{A2}
\]
For example the Bell-polynomial bounds are first applied with the \(E_r\), then with the saved Gaussian derivative bounds. Each term contains at least one derivative of \(C^k\), giving \(\lambda^k\); terms with several such factors are smaller.

Consequently the varying part of the central derivative satisfies the valid global refinement
\[
|D_K^{(r)}|\le18N^3 A_r
\sum_{k<K}\lambda^k
\le{18N^3 A_r\over1-\lambda},\qquad r=1,\ldots,4. \tag{A3}
\]
The zeroth-order cap bound \(|D_K|\le36KN^2\) is separate and is not removed by (A3). The new constants \(A_r\) can themselves be conservative; (A3) alone does not quantify a finite-sample panel reduction.

## 3. A full-panel bound using observed polygon slopes

**PROVED HERE.** Fix a closed panel \(I\) after removing all knot enclosures and splitting the weight jump at \(-7/4\). Suppose every varying probability query \(p_{\ell k}(t)\) stays in one known polygon segment throughout \(I\). This must be established by an enclosure of its whole image, not by inspecting its midpoint.

On that segment,
\[
Q_a(p)=a_{ai}+s_{ai}p,\qquad
s_{ai}=n(Y_{a(i+1)}-Y_{a(i)})\ge0.
\]
For the constant polygon part below \(1/n\), take \(s_{ai}=0\).
The slopes are functions of observed sorted heights. They introduce no latent or oracle information.

Let \(p_{\ell k,r}(I)\) be a rigorous upper bound on
\(\sup_I|p_{\ell k}^{(r)}|\), obtained from known-function interval derivatives or from (A2). Let \(Z_0(I)\) enclose the absolute value of the complete signed \(D_K\), including its cached reference term. For \(r\ge1\), take
\[
Z_r(I)=\sum_{k<K}\sum_{\ell=1}^{18}
|s_{a(\ell),i(\ell,k,I)}|\,p_{\ell k,r}(I). \tag{A4}
\]
The reference term is constant and contributes no derivative. Cancellation may be retained in a rigorous interval evaluation for \(Z_0\); it must not be assumed from opposite signs alone.

For \(V_r(I)\ge\sup_I|P_{M,c}^{(r)}|\), a complete central bound is
\[
C_{4,\mathrm{cen}}(I)=
\sum_{r=0}^{4}{4\choose r}V_{4-r}(I)Z_r(I). \tag{A5}
\]
This includes the term \(V_4Z_0\), which would be missed by bounding slopes alone.

For the outer integrand \(UQ_0(\Phi z)+WQ_3(\Phi z)\), use the same product formula separately for each summand. The order-zero polygon bound is its observed height range on the panel; its \(r\)-th derivative bound for \(r\ge1\) is
\[
|s_{ai}|\,\sup_I|\phi^{(r-1)}|.
\]
Add the two resulting fourth-derivative bounds. These formulas are valid for the rounded polygons as well, with the rounding error treated separately by the numerical-stability result.

**STANDARD RESULT APPLIED.** A two-node Gauss rule on a smooth panel of length \(\ell\) has absolute truncation error at most
\[
C_4(I)\ell^5/4320.
\]
Let \(S\) be a public upper bound on total integration length and allocate total Gauss error \(\epsilon_{\rm quad}\). Accept only panels satisfying
\[
C_4(I)\ell^5/4320\le
\epsilon_{\rm quad}\ell/S. \tag{A6}
\]
Adding (A6) over disjoint panels proves the total error budget. Subdivision with newly recomputed whole-panel bounds is allowed; an ordinary adaptive quadrature difference is not a replacement for (A5) or (A6). The original global envelope remains a fallback ensuring termination and the prior work upper bound.

Knot-enclosure gaps retain their own absolute-integral allocation. Their bound can also be localized: for each component \(G\) of the union, use a rigorous integrand bound \(B_0(G)\) and require
\(\sum_G B_0(G)|G|\le\epsilon_{\rm gap}\).
Overlapping enclosures are merged before this allocation. Node-location, primitive and summation errors still consume their separate shares.

## 4. What this changes about the saved panel counts

**PROVED HERE.** The saved counts can be reduced in a certificate whenever the corresponding observed/known local bounds are smaller:

- Replace \(N^3\) by the actual segment slopes in (A4); typical slopes are not silently assumed.
- Use the decaying known-map derivatives in (A2), rather than charging every central level equally for derivatives.
- Use whole-panel known-function bounds in (A5), and a rigorous signed bound for the actual central value, instead of global ranges covering all lattice arguments and coefficients.
- Merge overlapping knot gaps and split only at distinct enclosure components. The existing public knot count is already linear in \(N\); no smaller count is guaranteed merely by renaming the knot construction.

**OPEN.** No sampled polygons exist in this pass, so their slopes and resulting certified panel counts have not been measured. Extreme capped datasets can attain the same order as the global slope cap. Large finite constants in Bell-polynomial envelopes, knot localization, and central evaluation counts remain. The saved millions-of-panels figures are sufficient overcounts, not lower bounds; the refinements above neither establish an affordable runtime nor license substituting a guessed panel count.

**NUMERICALLY CHECKED, historical record only.** The previously saved public-count diagnostic reported 5,821,947 panels at \(N=1,000\) and 11,949,733,321 at \(N=1,000,000\), before the possible dyadic-mesh factor. Those records were read without rerunning them. They used global envelopes and are not measurements of a full estimator.

No complete interval backend was implemented to obtain another asymptotic claim. Whether implementation is scientifically justified depends on the complete influence-scale results and the [feasibility decision](feasibility_decision.md), separately from the static existence proof.
