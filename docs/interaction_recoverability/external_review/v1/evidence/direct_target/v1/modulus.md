**OPEN — sharp complete-experiment target modulus.** For fixed positive delta and all positive allocation proportions, define

`omega_(delta,pi)(epsilon)=sup{|I(eta)-I(eta')|: eta,eta' in the original class, sum_a pi_a H(P_eta,a,P_eta',a)^2<=epsilon^2}`,

with `H^2=integral(sqrt(p)-sqrt(q))^2`. Both profiles and all four coefficients vary in this supremum. This page proves a lower bound and a nonmatching constructive upper bound. It neither infers Hellinger closeness from DKW nor claims that a two-point lower bound provides an upper estimator.

**STANDARD RESULT APPLIED / PROVED HERE — exponential subclass lower bound.** Let `m=sqrt(e)`, `D=e-m`, `kappa=exp(.5)D`, `C0=49.3625+20.38 sqrt(2/pi) (approximately 65.6233873491624)`. In the admissible subclass g1=g2=exp, take amplitude coefficients `(u1,v1)=(m,e)` and `(u2,v2)=(b,m+e-b)`, with centered b alternatives in `[m,e]`. These are amplitudes, so the corresponding original shifts are their logarithms in `[.5,1]`. The audited score calculation gives

`d_pi(eta_b,eta_b') <= |b-b'| delta sqrt(pi2 C0)/(2m)`,

`|I_b-I_b'| >= kappa |b-b'|`.

Only source e2 changes; the remaining four laws have zero distance, not missing likelihood factors. Reverse the roles for e1. Hence, with `pi_dagger=min(pi1,pi2)`,

`omega_(delta,pi)(epsilon) >= kappa min{D, 2m epsilon/[delta sqrt(pi_dagger C0)]}`.

This pair stays in the original class and proves a global lower bound for the whole unknown-profile model. At these fixed boundary subclass centers it also gives a one-dimensional local lower direction; it is not a claim about every interior center. The known-profile theorem supplies its own matching subclass upper rate, which does not transfer to unknown profiles.

For finite counts use `pi_a=n_a/N`. Affinity multiplication gives `H_product^2<=sum_a n_a H_a^2=N d_pi^2`; also `TV_product<=H_product`. Choosing `d_pi<=1/(2 sqrt(N))` and thresholding any target estimator yields a constant times squared target separation as a risk lower bound. Thus the retained order is `c min{1,1/(min(n1,n2)delta^2)}`. If n1 or n2 is zero, that stratum's alternative is exactly unobserved and the lower risk remains constant. These statements concern all observations and use no Wasserstein-to-testing implication.

**PROVED HERE — valid Hellinger bridges.** For each a,

`||C_eta,a-C_eta',a||_infinity <= TV_a <= H_a <= epsilon/sqrt(pi_a)`.

Also, if both positive response laws have second moment at most M2bar_a, then

`W1(P,Q) <= integral y |p(y)-q(y)|dy <= 2 sqrt(M2bar_a) H(P,Q)`.

The first inequality follows from the dual Lipschitz representation, subtracting its value at zero. The second follows by writing `p-q=(sqrt(p)-sqrt(q))(sqrt(p)+sqrt(q))` and applying Cauchy–Schwarz; the other squared integral is at most four times the second-moment ceiling. Therefore, for a fixed signed bounded probability weight w,

`|integral w(p)[Q_P(p)-Q_Q(p)]dp| <= ||w||_infinity W1(P,Q)`.

This bounds direct weighted-quantile functionals from Hellinger distance in the required direction. It is not a bound of Hellinger distance by Wasserstein or CDF distance.

**PROVED HERE — explicit, nonsharp global upper bound.** Put `pi_*=min pi_a`, `t=epsilon/sqrt(pi_*)`, and choose coefficient truncation K_c, target truncation K and taper B. Let R,q,gamma be the coefficient-separation constants in estimator.md, `T=R+K_c+2`, `U=T+1`. For t small enough that `t<=Phi(-T)-Phi(-T-1)`, define

`E=t * 2e[exp(U)+exp(h(U))h'(U)]/phi(U)`,

`rho_delta(t,K_c)=min{.5, [(4K_c exp(R)+exp(R))E+4e(1+q)exp(-K_c)]/gamma}`.

To justify this comparison bound without an oracle coefficient assumption, compare the anchor telescopes of eta and eta'. Their quantile differences are at most E throughout the query range; thus profile differences at each required argument x are at most `2K_c E+4exp(x-K_c)`. In the two scaled coefficient rows, the sum of the two profile errors is at most `4K_c exp(R)E+4e(1+q)exp(-K_c)`. The observed e_i response difference adds `exp(R)E`. Apply P2 to coefficients c_i and c_i' with the first model's profiles held fixed. This proves `||c-c'||_infinity<=rho_delta`, with all five laws accounted for.

Now compare the two direct functionals using the *same* coefficient vector c in their weights. The private-anchor response laws under eta' determine its profiles independently of c'. The coefficient adjustment is at most C_c rho_delta by the global target Lipschitz bound. Both models have a telescope and taper error; the direct quantile comparison uses the preceding Hellinger–Wasserstein bound. Hence

`omega_(delta,pi)(epsilon) <= C_c rho_delta(t,K_c) + 2 Rplus(K) + 2 Tplus(B) + 2 epsilon [sum_(a=0,3,4) Mplus(K,B)^2 M2bar_a/pi_a]^(1/2)`.

The right side can always be capped by 2Mbar. All constants are defined and bounded explicitly in estimator.md; none is a fitted nuisance radius. This is a direct-functional upper bound with separate coefficient error, rather than target-stage reconstruction of whole profiles.

For fixed delta and pi, let `L=log(sqrt(pi_*)/epsilon)`, `K_c=K=max(1,floor(.5 sqrt(L)))`, and `B=max(1,.5 sqrt(L))`. As epsilon tends to zero the quantile condition eventually holds. The coefficient observation term is `exp[-(1-(.5+delta)/4)L+o(L)]`; the direct weighted-quantile term is `exp[-7L/8+o(L)]`. Their decay is faster than the telescope, coefficient-truncation and taper terms. Consequently

`omega_(delta,pi)(epsilon)=O_(delta,pi)(exp[-.5 sqrt(log(sqrt(pi_*)/epsilon))])`.

This converges to zero uniformly over the original parameter class at fixed delta and pi, but is much larger than the linear lower order. R grows as delta^(-1/2); no fixed-delta threshold or constant is silently asserted uniform in delta. The existing root-(n delta^2) subclass lower order remains the only sharp separation dependence established here.

**OPEN — the single sharp inequality this attempt did not prove or refute.** Does there exist, for every fixed `delta>0` and `pi_a>0`, a finite constant C_(delta,pi) and epsilon0 such that

`|I(eta)-I(eta')| <= C_(delta,pi) d_pi(eta,eta')`

for **all** original-class eta,eta' with `d_pi<=epsilon0`? The lower bound makes linear order the first candidate to check, but does not establish that it is achievable. The direct upper bound above does not resolve this inequality. If it is false, a stronger admissible complete-five-law lower construction and an appropriate slower attainable upper order would still be needed. No such construction was obtained in this bounded attempt.

**STANDARD RESULT APPLIED / OPEN — local versus global.** The global supremum permits profiles and coefficient centers to drift with epsilon. Its behavior need not equal the local modulus at a fixed parameter. At a fixed interior parameter, the usual score criterion asks whether the target derivative has a finite-norm representer in the complete weighted score space. For a sufficiently regular candidate influence, let `k_a(z)=pi_a psi_a'(F_a(z))phi(z)`, and set `A_a=(0,alpha1,alpha2,1,0)`, `B_a=(0,beta1,beta2,0,1)`. A sufficient profile adjoint system is

`sum_a k_a(x-A_a)=K1(x)`,

`sum_a k_a(h_inverse(x-B_a))/h'(h_inverse(x-B_a))=K2(x)`.

It must also satisfy all four coefficient moments

`integral k_i(z)g1'(z+alpha_i)dz = partial_(alpha_i) I`,

`integral k_i(z)g2'(h(z)+beta_i)dz = partial_(beta_i) I`, `i=1,2`,

and centered primitives `integral^z k_a(s)F_a'(s)/(pi_a phi(s))ds` must lie in L2(phi). More generally formulate these equations distributionally, allowing annihilators of the normalization constraint `r_j(0)=0`; restricting to smooth k without this qualification could wrongly exclude nonsmooth influence functions. The general criterion is the equality against every admissible tangent, as in historical P4.

The anchor-only identity sets k1=k2=0 and therefore misses generic nonzero coefficient derivatives. Its periodic-tail infinite variance, proved in target_identity.md, is an obstruction for that representation even with known coefficients. It neither proves the all-five-law inequality false nor rules out a finite-norm local representer using e1/e2. Neither that representer nor a countersequence disproving it was constructed. Even a local solution would still require nonlinear and estimation-error control for an attainable estimator; it would not automatically establish the global inequality.

**STANDARD RESULT APPLIED / OPEN — confidence sets and likelihood.** Five DKW CDF bands cover the true source laws at their allocated failure level. An original-class target interval formed by the infimum and supremum of I over that set would have coverage, but its radius is governed by a CDF-compatible target modulus. DKW supplies no Hellinger radius, and the present Hellinger bound cannot be inserted in that direction. A target profile-likelihood or Hellinger testing construction can use the same five laws and profile out nuisance parameters, but requires entropy/testing/coverage and approximation/optimization proofs in this infinite-dimensional class. A dimension-free root-N Hellinger ball is not established. Finite-sieve feasible endpoints are inner witnesses and cannot certify full-class extrema. The direct estimator's sampling upper bound is instead proved directly for its integrated empirical quantiles, so it does not rely on either missing bridge.
