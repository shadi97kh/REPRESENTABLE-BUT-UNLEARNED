**PROVED HERE — exact finite direct-target identity.** Under [the original model](exact_problem.md), let `f1=phi` and `f2(x)=phi(h_inverse(x))/h'(h_inverse(x))`. For component j with shifts `(c1,c2)=(alpha1,alpha2)` or `(beta1,beta2)`, define

`K_j(x;c)=f_j(x-c1-c2)-f_j(x-c1)-f_j(x-c2)+f_j(x)`.

For every admissible pair of profiles, every coefficient vector in the original box and every integer `K>=1`,

`I12 = sum_j integral K_j(x) g_j(x) dx`,

`D1(x)=q3(x)-q0(x)`, `D2(x)=q4(h_inverse(x))-q0(h_inverse(x))`,

`W_(j,K)(x)=sum_(k=1)^K K_j(x+k)`,

`I12 = sum_j integral W_(j,K)(x) D_j(x) dx + R_K`,

`R_K=sum_j integral K_j(x) g_j(x-K) dx`.

Proof: change variables from the component index in each of the four target terms to the profile argument. In particular `x=h(z)` contributes the Jacobian `1/h'(h_inverse(x))`; the density is not `phi(h_inverse(x))` alone. The private anchors give `D_j(x)=g_j(x+1)-g_j(x)`. The finite telescope is `g_j(x)=sum_(k=1)^K D_j(x-k)+g_j(x-K)`. Substitute into the first display and set `u=x-k` in each of the finitely many integrals. This produces `K_j(u+k)`, fixing both the direction and sign of the weight shift. All terms are integrable by the exponential profile envelopes and the known exponential moments of f_j.

**PROVED HERE — cancellation and remainder.** Both `integral K_j=0` and `integral x K_j(x) dx=0`; the target annihilates affine profile perturbations. Before taking absolute values, rewrite its remainder as

`R_(j,K)=E integral_0^c1 integral_0^c2 g_j''(X_j-K+s+t) ds dt`,

where `X1=Z`, `X2=h(Z)`. Consequently

`|R_K| <= 10 exp(-K) sum_j m_(1,j)(exp(c1)-1)(exp(c2)-1)`.

One can take the minimum with `2 exp(-K) sum_j m_(1,j)(1+exp(c1))(1+exp(c2))`, obtained directly from the profile envelope. The derivative bound retains signed second-difference cancellation, but it does not always improve the constant throughout the coefficient box. Neither bound establishes a new rate.

**PROVED HERE — justified infinite identity.** Since `0<D_j(x)<=2(e-1)exp(x)`,

`sum_(k>=1) integral exp(x)|K_j(x+k)| dx = [1/(e-1)] integral exp(u)|K_j(u)|du`

and the final integral is at most `m_(1,j)(1+exp(c1))(1+exp(c2))`. Thus the integral of the absolute series multiplied by `D_j` is finite. Tonelli followed by dominated convergence justifies interchanging the infinite sum and integral in the direct identity. Also `R_K -> 0` uniformly in the original class at rate at most a constant times `exp(-K)`. These facts establish existence of the infinite functional, not finite variance of its empirical-quantile estimate.

**PROVED HERE — shared reference weight.** Returning to the latent quantile coordinate gives

`A3(z)=W_(1,K)(z)`, `A4(z)=h'(z) W_(2,K)(h(z))`, `A0(z)=-A3(z)-A4(z)`.

The truncated target `I12-R_K` equals `sum_(a=0,3,4) integral A_a(z) q_a(z) dz`. In particular the two occurrences of q0 are the same estimated quantile curve. They must be combined before evaluating a variance. There is no independence assumption between those two contributions. The e1/e2 observations enter through coefficient learning; treating the displayed known-coefficient identity as an estimator for unknown coefficients would omit their contribution.

**PROVED HERE — taper bound uniform in K.** Let `chi_B=1` for `|z|<=B`, zero for `|z|>=B+1`, and `chi_B=1-3t^2+2t^3` in between, with `t=|z|-B`. This is C1 and `|chi_B'|<=1.5`. Multiplying the direct functional by chi_B incurs error at most

`T_(K,B)=2(e-1) integral (1-chi_B(z)) [|W_(1,K)(z)|exp(z)+h'(z)|W_(2,K)(h(z))|exp(h(z))] dz`.

The following fully known uniform bound is sufficient:

`T_(K,B) <= C_T(delta) exp(-B)`,

`C_T(delta)=2(e-1)[16+((1+e^2)^2/(e^2-1))(m_(2,1)+m_(2,2))]`.

Proof: both f_j are symmetric unimodal densities. For f2, phi decreases and h' increases on the positive half-line. The sum of a density over a unit lattice is at most `1+TV(f_j)=1+2f_j(0)<2`, by comparing each interval integral with its endpoint value and summing the variation. Each weight is a sum of four shifted lattice sums, so `|W_(j,K)|<=8`. The negative x-tail integral of `exp(x)|W|` is therefore at most `8 exp(-B)`, and for component 2 the cutoff is `h(-B)<=-B`. For the positive tail, one shifted term obeys

`integral_(x>L) exp(x) f_j(x+k-tau) dx <= m_(2,j) exp(-L+2tau-2k)`.

Sum k and the four shifts tau to obtain the coefficient `(1+e^2)^2/(e^2-1)`. The component-2 positive cutoff is `h(B)>=B`. Combining the two negative and two positive tails proves the display. This controls the directly integrated functional without bounding reconstructed profiles over extreme quantiles.

**PROVED HERE — failure of the naive untruncated influence function, restricted scope.** At the admissible exponential profiles and `alpha1=alpha2=1/2`, Gaussian periodization yields

`W_(1,infinity)(-m+r) -> P(r)=8 sum_(k>=1, k odd) exp(-2 pi^2 k^2) cos(2 pi k r)`

as integer m tends to infinity. The Fourier expansion follows by periodizing phi (whose k-th Fourier coefficient is `exp(-2 pi^2 k^2)`) and taking its translated second difference. The first coefficient is nonzero. More generally the coefficient at integer k is multiplied by `(exp(-2 pi i k alpha1)-1)(exp(-2 pi i k alpha2)-1)`, so it is nonzero at k=1 whenever both alpha shifts are strictly below one.

For the particular infinite anchor-only L-statistic, its source-e3 influence primitive has derivative in latent coordinate

`d ell3(F3(z))/dz = W_(1,infinity)(z) F3'(z)/phi(z)`.

At exponential profiles `F3'(z)~exp(z+1)` as z tends to negative infinity. Endpoint integration on infinitely many fixed-length phase intervals on which P stays away from zero gives influence magnitude bounded below by a constant times `exp(z^2/2+z)/|z|`. To see endpoint dominance, integrate over the last interval of width order `1/|z|`; P changes by O(1/|z|), whereas contributions a fixed phase interval away are exponentially smaller. Integration by parts gives the same leading endpoint term. Its square against phi has divergent integral, of order `exp(z^2/2+2z)/z^2` on these intervals. Therefore this particular untruncated representation has infinite influence variance. Reference covariance cannot cancel independent e3 variance.

This does not prove all-five-law nonregularity or a global minimax obstruction. Alternative e1/e2 contributions and coefficient corrections can change the representer. If either alpha shift equals one, the periodic term vanishes; for example at both shifts one, `W_(1,infinity)(x)=phi(x-1)-phi(x)`. The counterexample diagnoses an unjustified step in the naive formula, not impossibility of the research objective.

**NUMERICALLY CHECKED.** [checks.json](checks.json) records two analytic fixtures, finite identities, source-coordinate Jacobians/reference combination, density normalization, inverse round trips, a primitive derivative, the periodic residue and the integer-shift exception. Largest identity error: `1.2795320358804929e-14`; source-weight error: `1.4210854715202004e-14`. The periodic residue at phase .125 was `1.513371417176337e-8`, versus Fourier value `1.5133714240924264e-8`. These float64 checks on stated finite integration intervals are not certified tail or roundoff enclosures.
