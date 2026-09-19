# Consequences, exact overlap and significance

The mathematical decision is [worst-case exponential efficiency established](decision.md), on the original profile class and the stated small-\(\delta\) warp interval. This is a model-specific analytic result, not a numerical discovery, a general ML principle or a proof of publication priority.

## 1. Which quantities have been resolved

**PROVED HERE.** The new construction gives an admissible truth \(\eta_{\delta,\theta}(0)\) for every fixed known \(\theta\in[1/2,2]\), \(0<\delta\le1/3{,}200{,}000\), with
\[
V_{\rm eff}(\eta_{\delta,\theta}(0);\delta,\theta)
\ge \frac{e^{11/5}}{32000}e^{1/(640\delta)}.
\]
Together with the saved uniform constructed-influence upper bound,
\[
\log\!\left(1+\sup_{\eta\in\Theta}
V_{\rm eff}(\eta;\delta,\theta)\right)=\Theta(1/\delta),
\]
uniformly in \(\theta\) as \(\delta\downarrow0\). This closes the previous question of whether exponential **worst-case local variance** could be solely an artifact of that particular estimator. A polynomial uniform bound for compatible regular influences on the original class is impossible.

The finite-experiment construction separately gives, for every \(N=5n\), \(n\ge1\),
\[
\mathcal R_{N,\delta,\theta}\ge
\frac{e^{11/5}}{2{,}048{,}000}
\min\left\{1,\frac{e^{1/(640\delta)}}N,
\frac1{\delta\sqrt N}\right\}.
\]
The constants and full nonlinear path argument are in [lower_bound.md](lower_bound.md).

| Quantity | Consequence of the new proof |
|---|---|
| Worst-case efficient variance over the original class | Exponential in inverse separation at the level of logarithmic order |
| Variance at the explicit hard truth | Explicit exponential lower bound; exact value unknown |
| Efficient variance at the historical \(\delta=.1\) exponential reference | Unchanged and unresolved |
| Variance of the previously constructed pointwise estimator | Its upper bound remains valid; efficiency is not established |
| Finite-sample minimax squared risk | Explicit three-regime lower bound; matching upper bound unknown |
| Sample size sufficient or necessary for a fixed practical tolerance | No exponential requirement follows just from the local variance constant |

The new lower bound is achieved by choosing profiles as a function of the known warp. It makes no assertion that every fixed profile pair has exponential local variance. In particular the old reference has both different profiles and a warp outside the new small-\(\delta\) interval.

## 2. Why the nonlinear term matters

**PROVED HERE.** The exact coefficient balance makes the score at \(t=0\) exponentially small, but away from zero the shifted profile derivative differs by a term of order \(\delta|t|\). Integrating the actual density score gives
\[
H(P_{1,t},P_{1,0})
\le20e^{-1/(1280\delta)}|t|+10\delta t^2.
\]
This second term cannot be discarded when taking testing separations depending on \(N,\delta\).
It produces the intermediate squared-risk scale \(1/(\delta\sqrt N)\) in the lower bound. Equating that scale with the local term gives
\[
N\asymp\delta^2e^{1/(320\delta)}
\]
up to constants for this proved inequality.

These are regimes of a lower bound, not a complete minimax characterization or demonstrated estimator performance. The fixed-\(\delta\) information regime may begin extremely late. Extrapolating its \(e^{1/(640\delta)}/N\) term to all sample sizes would overstate the proved lower bound.

A claim of exponential sample requirements for fixed accuracy would require an appropriate uniform finite-experiment argument at a fixed target separation. This path instead retains the quadratic likelihood movement. The result also supplies no polynomial sufficient sample bound: that would need an upper estimator analysis.

## 3. Trabs and the information reduction

**STANDARD RESULT APPLIED.** [Trabs, arXiv:1307.6610v2, Theorem 2.7](https://arxiv.org/html/1307.6610v2) supplies the adjoint-score regularity and convolution-information framework already used in the project. The present local lower bound is exactly the ordinary one-submodel inequality
\[
|DI[q]|^2\le\|\psi\|_\pi^2\|Sq\|_\pi^2.
\]
Only one admissible direction is required for a lower bound. This is not a finite approximation to the full tangent closure: every full-model compatible influence must satisfy that one equation as well.

The constructed direction has strictly positive information at every fixed \(\delta>0\), although the information is exponentially small. Calling it exactly singular at positive separation would be incorrect. The separate retained theorem supplies a finite compatible influence and rules out any inference of infinite variance at a fixed allowed warp.

The new content is the explicit admissible profile pair, exact four-law invariance, exponentially small score in the fifth law, and nonvanishing interaction derivative. The information inequality itself is established machinery.

## 4. Heinrich–Kahn: pointwise and uniform rates

**STANDARD RESULT APPLIED AS A COMPARISON.** The primary
[Heinrich–Kahn paper, arXiv:1507.04313v1](https://arxiv.org/html/1507.04313v1), treats finite mixtures of known component kernels with an unknown mixing distribution. Its loss is expected Wasserstein-\(1\) error for that mixing distribution.

- Theorem 3.2, under Assumption A, proves a local minimax lower rate \(n^{-1/[4(m-m_0)+2]}\).
- Theorem 3.3, under \(B(2m)\), provides the matching local upper rate and global rate \(n^{-1/(4m-2)}\).
- Theorem 3.5, under \(B(1)\), supplies a nonuniform pointwise \(O(n^{-1/2})\) result at every fixed finite mixing distribution.
- Section 3.3 explains the distinction through truth-dependent constants and delayed pointwise asymptotic regimes.

These precise theorem statements support the conceptual comparison. They do not supply any rate for the present model.

| Comparison | Heinrich–Kahn | Present construction |
|---|---|---|
| Observation model | Finite mixture of known kernels | Five monotone response pushforwards, with two unknown analytic profiles and known shared-latent warp |
| Unknown target | Mixing distribution | Scalar interaction contrast |
| Loss | Expected \(W_1\) distance | Ordinary squared error |
| Hardness mechanism | Nearby mixing distributions and higher-order mixture degeneracy | A saturated profile derivative closely matches a warped exponential, with an exact coefficient cancellation path |
| Local versus uniform issue | Pointwise parametric rate can coexist with slower local/global minimax rates | Fixed-warp regular attainment coexists with exponentially large worst-case information constants and a nonlinear testing regime |

The adjacent results show that weak/singular information and the distinction between pointwise and uniform recovery are established themes. Their mixture assumptions and exponents cannot be transplanted. The present result requires its own global profile admissibility and full nonlinear density calculations, which are supplied separately.

## 5. Exact additional contribution and limits

**PROVED HERE as a model-specific obstruction.** Smooth saturation limits the profile's derivative ratio globally while preserving an exponentially accurate cancellation on a Gaussian region whose radius grows as \(\delta^{-1/2}\). The tail probability controls the residual information. The profile is kept fixed along an exact nonlinear coefficient path, so four observed laws are identical and the fifth remains quantitatively close for an entire testing segment. The interaction nevertheless changes at a uniformly nonzero rate.

Combined with the preserved family attainment result, this establishes an infinite-dimensional example with:

1. Pointwise regular interaction recovery from exactly five source laws with both profiles and all four coefficients unknown.
2. An unavoidable exponential worst-case local information cost in inverse warp separation within that same original class.
3. A separate finite-experiment lower bound retaining the nonlinear intermediate regime.

Smooth saturation, density-score calculations, Cauchy–Schwarz information bounds, Hellinger product affinity and two-point testing remain established tools. The retained contraction and quantile estimator is also based on established inverse-functional machinery. Equally informed direct-target and inverse methods are subject to the same model lower bound, but their relative performance is not determined.

**OPEN:** publication priority and whether this model-specific obstruction plus attainment is sufficiently substantial for a paper. This bounded source comparison does not establish an absence of antecedents. It strengthens the significance case beyond merely reporting a poor constructed upper constant; it does not establish a new general learning principle.

The principal remaining statistical gap is a matching finite-sample minimax upper bound, including the separation dependence and nonlinear regime. Leading exponential constants, efficient influences/estimators and practical accuracy remain unresolved. The old reference-efficiency question, failed learned-target pilot, backend no-go and lack of biological siRNA validation remain unchanged.

## Source inspection

Read the retained warp-family comparison and upper-bound documents and the historical efficiency argument's score/norm definitions and fixed reference. Re-read cached Trabs v2, Theorem 2.7, in
runs/interaction_recoverability/prep-20260912-082154/unknown-profiles-v1/trabs_v2.txt, lines 330–390.
No cached Heinrich–Kahn text was found; its primary HTML was accessible and Sections 2.4, 3.1–3.3 and the adjoining theorem statements were inspected. No requested primary source was unavailable.

The paper comparison was performed by a separate agent inside the same system and integrated into this static self-review. It is not external independent peer review. No numerical experiments or computational allowance operations were performed.
