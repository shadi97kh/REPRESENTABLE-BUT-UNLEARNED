# Algorithm, operation budget and implementation boundary

**PROVED HERE:** the mathematical accelerated statistic is specified by the following steps. All arguments are observed responses or public constants; true profiles, coefficients, Jacobians, influence values and latent samples are absent. [Theorem](theorem.md) proves its statistical expansion and [quadrature proof](quadrature_proof.md) gives the numerical allocation.

1. Check five balanced groups, \(N=5n\ge10\), finite positive responses and the observable cap \(Y\le N^2\). The mathematical rule returns zero on failure. On the original cap event, upward-round responses to the specified dyadic input accuracy and clip rounded values at \(N^2\). Sort within each source and create polygons at \(i/n\).
2. Compute public \(K,R,M\), centers, rectangles and preconditioners. Construct \(E_{45}\), \(C\), \(D_K\) and the normalized profile evaluation operators. Only coefficient calibration uses the finite transported profile values.
3. Starting at each public center, take \(2\lceil\log N/|\log(3/4)|\rceil\) projected updates. Use the correct \(-1,-3\) calibration charts, whose entire domains have been verified. Do not select a true-root branch.
4. Form the finite signed target weights and the finite central difference integral. Cache the reference-point part of \(D_K\).
5. Enumerate all polygon-knot preimages for each varying central probability map and the outer \(\Phi(z)\) map. Use forward certified bisection. Merge overlapping enclosures, exclude the union with its absolute error bound, and split at \(b=-7/4\).
6. Integrate each remaining smooth interval by two-node Gauss panels using the **proved** fourth-derivative envelope, including observed polygon slopes. Apply the allocated endpoint, primitive and summation tolerances. Clip the final target at 1000.

These steps are a fully specified certifiable mathematical algorithm. The delivered `estimator.py` supplies finite anchor/profile/weight operators, public projected iteration, diagnostic monotone bisection, a two-node Gauss component with an explicit derivative-bound argument, and an unconditional sampled-data refusal. It is not a full certified interval implementation or an end-to-end execution. It intentionally does not pretend that an adaptive quadrature error estimate is a certificate.

## Explicit sensitivities and allocation

The companion proof defines \({\cal C}_4,B_0,J_c,J_Y\) through finite positive Bell-polynomial formulas. Here are usable conservative coefficients for the iteration/input part of that graph.

Every calibration profile query has derivative at most \((36K+7)N^3\) on its fixed chart. A row of \(\widehat B\) uses at most two such terms. Thus
\[
C_B=100,\quad \operatorname{Lip}(\widehat B)\le100(K+1)N^3.
\]
The public matrices have \(\|M_{i0}^{-1}\|\le69\): their Frobenius condition product is below 69 and \(\|M_{i0}\|_F>1\).
Consequently one may take
\[
L_N=1+6900N^3(K+1),\qquad
A_T=\sum_{s=0}^{T-1}L_N^s,\qquad C_{\rm data}=10350(K+1).
\]
For the latter, a response-height perturbation of size \(\tau\) changes one residual row by at most \((72K+18)\tau\); bound the vector norm by twice that amount and multiply by 69.

Using the notation \(H_1,{\cal R}_1,P_0,U_0,W_0\) of the quadrature proof and \(B=R+M+4\), set
\[
L_c=2+2H_1(B){\cal R}_1,
\]
\[
J_c=2N^2\{2R[2H_1(R){\cal R}_1+2ML_c]+36K(2M+1)L_c\},
\quad J_Y=2R(U_0+W_0)+36KP_0.
\]
Each coefficient derivative of \(K_j\) is a difference of two first density derivatives; these formulas therefore bound the complete four-vector target sensitivity and the fixed-coefficient height sensitivity. They need no unknown derivatives.

Give **each pair** final coefficient error at most
\(\delta_c=[16N(1+J_c)]^{-1}\). The complete four-vector error is at most \(2\delta_c\). Bound every computed update error by \(\delta_c/A_T\), allocating one quarter to input rounding and the rest to primitives, preconditioner and arithmetic. Set
\[
\tau_Y\le\min\{\delta_c/(4A_TC_{\rm data}),[32N(1+J_Y)]^{-1}\}.
\]
Gap exclusion, Gauss truncation, integration node/weight errors, and final summation/endpoint errors each receive \(1/(8N)\). The coefficient target share is at most \(1/(8N)\), and the fixed-coefficient input share at most \(1/(32N)\). Together these are below \(1/N\); remaining slack can cover outward rounding of known constants. All graph enclosures must meet the allocated shares. A uniform floating-point precision setting alone is insufficient.

The conservative comparison permits worst-case empirical iteration amplification. It requires \(O((\log N)^2)\) fractional bits, not a fixed \(N^{-D}\) input tolerance. A separate stochastic stability argument might improve this sufficient precision, but no such improvement is used in the complexity or count diagnostics.

**Schedule comparisons in finite arithmetic.** Select \(K\) and \(k_N\) by exact rational-power comparisons \((271/290)^j\le N^{-2}\) and \((3/4)^j\le N^{-1}\), respectively. This avoids undecidable exact logarithm-ceiling comparisons. Enclose \(R\) to its allocated endpoint error; choose the lattice integer from the upper rational endpoint. It differs from the displayed ideal \(\lceil4R+30\rceil+2\) by at most one once the enclosure width is below \(1/4\). The added lattice shell has the already proved \(N^{-10+o(1)}\) error on capped data; it cannot change the root-\(N\) comparison. Choose a dyadic mesh below the proved Gauss limit using rational upper derivative bounds. Polygon knot decisions at equality use enclosures, never an exact comparison of transcendental quantities.

## Work and memory

Let \({\cal B}\le36n/(1-\lambda)+36K+2n+4\) bound knot preimages and \(S=1+2R\). Then
\[
h=\min\{1,[4320/(8N^4{\cal C}_4S)]^{1/4}\},\qquad
J_{\rm panels}\le{\cal B}+\lceil S/h\rceil+2.
\]
A dyadic mesh within a factor two increases only the corresponding mesh-panel bound by at most two.

| Stage | Mathematical work bound | Information used |
|---|---|---|
| Input and sorting | \(O(N\log N)\) comparisons | Observed source heights at specified increasing precision |
| Coefficient iteration | \(O(\log^2N)\) quantile/inverse evaluations | Five source polygons and public preconditioner |
| Preimage localization | \(O({\cal B}K\log N)\) primitive evaluations | Known monotone maps and knots \(i/n\) |
| Target panels | \(O(J_{\rm panels}(K+M))\) quantile/inverse evaluations | Finite signed weights, polygons and reference cache |
| Storage | \(O(N+K+M)\) words with streamed panels | Sorted observations and knot enclosures |
| Arithmetic | \(O((\log N)^2)\) bits; each primitive polynomial in this count | Certified range-reduced elementary functions |
| Total | \(N\,\mathrm{polylog}N\) bit operations, \(N\,\mathrm{polylog}N\) bits of memory | Includes reading the required input precision |

This improves the saved mathematical \(N^{6+o(1)}\) implementation bound. It does not prove optimality or a small multiplicative constant.

## Actual finite-count assessment — no estimator evaluated

**NUMERICALLY CHECKED:** [checks.json](checks.json) and [practical_counts.json](practical_counts.json) evaluate only fixed analytic functions and public count formulas. Total \(N\) is divided equally among five groups.

| Quantity | \(N=1{,}000\) | \(N=1{,}000{,}000\) |
|---|---:|---:|
| Central depth \(K\) | 206 | 410 |
| Updates per coefficient pair | 50 | 98 |
| Lattice \(M\) | 45 | 51 |
| Response radius \(R\) | 3.21895 | 4.55228 |
| Minimum central tail probability | \(5.77025\,10^{-4}\) | same |
| Corresponding source effective rank \(n p\) | about .115 | about 115.4 |
| Tail guard probability \(\Phi(-(R+1))\) | \(1.22722\,10^{-5}\) | \(1.40983\,10^{-8}\) |
| Tail guard source effective rank | about .00245 | about .00282 |
| Public knot overcount | 117,715 | 110,309,501 |
| Public Gauss panel overcount | 5,821,947 | 11,949,733,321 |
| Sufficient input fractional bits from allocated \(\tau_Y\) | 2,550 | 7,994 |

These panel counts use the real mesh formula, before the possible factor-two dyadic-mesh overhead. They are evaluations of conservative sufficient formulas, not lower bounds, measured runtimes or estimates of the work an optimized local interval method would actually need. The coefficient count is per pair. A central value has at most \(36K\) response evaluations before caching its reference part; the cache removes half on subsequent calls.

The central probability requirement has improved from the old \(\Phi(-7)\). The response-tail requirement is separate and remains poor at these example counts. Neither sufficient condition is an operational zero-output gate or a necessary sample size. In particular the old astronomical proof-onset diagnostic has not been used to choose an experiment.

**OPEN / practical objective not reached:** the conservative full certificate is expensive even at small \(N\), and the statistical constants do not demonstrate an informative finite-sample regime. A complete interval backend and sharper local integration/empirical-error bounds would be needed before describing an affordable certified validation. Ordinary adaptive quadrature could produce exploratory numbers, but would not establish the theorem's numerical obligations. This continuation therefore identifies no justified sampled-validation campaign and requests none.

