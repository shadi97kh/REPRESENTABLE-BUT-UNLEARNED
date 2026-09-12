> Publication copy of `docs/interaction_recoverability/learned_target_v1/theorem.md`. Mathematical and result content is retained; execution material is omitted. This dated record retains its original conclusions and limitations, including conclusions superseded by later stages. Original SHA-256: `d2d36cd66d01c5a9e6220438836e617f3a7e1e2bc4dfc544d893de4c551e0acd`.

# Proof audit: learned correction in the nonlinear experiment

The proposed sharp learned-estimator guarantee is **OPEN**. The completed results below are exact risk accounting, standard finite-dimensional regularization, and a coarse nonshrinking bound. None proves the desired learned nonlinear rate or resolves historical P6.

## 1. Fixed-outcome score and target derivative — PROVED HERE

For a differentiable admissible path eta_t, at fixed observed y write `z_t=F_t^{-1}(y)`. Implicit differentiation gives `dot z=-dot F/F'=-u`. Since `log p_t(y)=log phi(z_t)-log F_t'(z_t)`,

`dot log p(y)=z*u-dotF'/F'+F''*dotF/(F')^2=z*u-u'`.

Integration by parts against the Gaussian density gives mean zero when `phi(z)*u(z)` vanishes at both ends. The finite analytic direction bank has `|u|<=4||v||_1` by the derivative envelopes, so this boundary condition and score square integrability hold. Differentiating the four-action target under its integrable exponential envelope gives the implemented b. Finite differences at fixed y, density normalization, score normalization and independent target derivatives check these identities in the correctness suite.

## 2. Exact conditional risk — STANDARD RESULT APPLIED

Condition on the fit and selection data for **one rotation**, including the learned model, directions and lambda. Write J for the plug-in value actually used (it may itself be a numerical integral) and psi for the fixed correction actually evaluated. For arbitrary deterministic weights pi define

`B0 = I0 - J - sum_a pi_a E0,a[psi_a(Y)]`.

Independence of the evaluation observations gives exactly

`E[(Ihat-I0)^2 | F,C] = B0^2 + sum_a pi_a^2 Var0,a(psi_a)/n_a,E`.

The estimator's bias is **-B0**, so the correction has a plus sign. Taking `pi_a=n_a,E/N_E` reduces the variance to `sum_a pi_a Var0,a(psi_a)/N_E`. Exact fitted-law centering is needed to interpret B's value and derivative at the fitted center, but is not needed for this general MSE identity: numerical centering error belongs to B0. Fitted-law variance is not substituted for true-law variance without an additional transfer bound.

For the average of three dependent rotations, pointwise `(mean_k error_k)^2 <= mean_k error_k^2`. Expectations therefore give a valid, possibly loose mean-of-risks bound. There is no factor-three variance improvement claim. The recorded per-fold empirical variances are descriptive estimates, not calibrated cross-fit confidence intervals.

## 3. Oracle linearized calculation — STANDARD RESULT APPLIED

Fix a finite direction space with exact positive metric M and whiten it. In the local linear experiment, suppose a source-score statistic has mean `G theta` and covariance `G/N_E`, and the target increment is `b^T theta`. For a fixed lambda, the correction is `c=(G+lambda*Id)^(-1)b`. In a G-eigenbasis with eigenvalues kappa_j,

`Var = (1/N_E) sum_j kappa_j*b_j^2/(kappa_j+lambda)^2`.

Its residual target is `R=b-Gc`, with coordinates `lambda*b_j/(kappa_j+lambda)`. Cauchy–Schwarz, with equality in the residual direction, yields

`sup_{||theta||<=r} Bias^2 = r^2 sum_j lambda^2*b_j^2/(kappa_j+lambda)^2`.

These are exact statements about this linear experiment, for fixed lambda and known geometry. A zero kappa with nonzero b contributes its full squared target component to the bias. Ridge does not identify that component. The correctness fixture includes this case. Replacing the oracle r by `.02`, fitting G, or choosing lambda from a proxy does not preserve an oracle inequality automatically.

The finite adaptive search is also exact standard algebra: on `v^T M v=1`, the objective is `(R^T v)^2 / (v^T(G+tau*M)v)`. Its maximum is `R^T(G+tau*M)^(-1)R` and a maximizing vector is proportional to `(G+tau*M)^(-1)R`. This is a generalized Rayleigh quotient followed by Galerkin enrichment. In the actual implementation M is numerical; the exact algebra applies to that matrix. It establishes neither an infinite-class maximum nor a new learning principle.

## 4. Finite nonlinear paths — PROVED HERE / NUMERICALLY CHECKED

Let eta_t be the additive profile/shift path, with finite t restricted by the globally proved envelope and coefficient margins in `method.md`. With psi fixed and exactly fitted-law centered, differentiation of the whole source law gives

`B(0)=0`, `B'(0)=DI[v]-sum_a pi_a E_hat[psi_a*s_a,v]=R(v)`.

Where the second derivative exists and is integrable, the exact remainder formula is

`B(t)=t*R(v)+integral_0^t (t-u)*B''(u) du`.

The signed integral convention handles negative t. This formula includes coefficient/profile coupling through changed indices. The numerical fixture confirms the linear derivative and detects a nonzero second-order remainder. Actual fitted-model paths store B, its linear prediction, their difference and source compatibility. Neither a small observed remainder nor failure to find a large residual is a uniform upper bound on the original parameter class. Finite quadrature and centering introduce an additional numerical baseline error in the code's path comparison.

## 5. A uniform but uninformative learned-correction bound — PROVED HERE

The following demonstrates finite risk under the actual global envelopes without assuming full-profile consistency. It does **not** establish a shrinking rate. Work in real arithmetic and condition on any finite fitted coefficient vector c, finite quadrature plug-in J, and fitted-law centering rule supported on `[-L,L]` with nonnegative weights summing to one. Here L=10. Let V be the L1 norm of c in the physical 36-direction bank.

The direction atoms obey `|k|<=2`, `|k'|<=3`. The source coefficient directions obey `|g'|<=2 exp(x)`, `|g''|<=10 exp(x)`. Hence `|dotF/F'|<=4V` and `|dotF'/F'|<=20V`. Also `1<=h'<=1.1+.2|z|`, `|h''|<=.2`, and

`|F''/F'| <= 20*h' + .2 <= 22.2+4|z|`.

Thus `|s_c(F_hat(z))| <= 108.8V+20V|z| <=109V(1+|z|)`.

For any original-class truth and fitted original-class model, at the same source action and latent z the two response values have ratio between `1/(4*exp(.5))` and `4*exp(.5)`. The `.5` is the maximum single-action shift difference; private anchor shifts are fixed. Moreover `(log F_hat)' >=1/4`. Therefore

`|F_hat^{-1}(F_0(z))-z| <= D := 4*log(4*exp(.5))`.

Numerical centering supported on `[-L,L]` is bounded by `109V(1+L)`. With `A=2+L+D` and `m=E|Z|=sqrt(2/pi)`,

`E0 psi_c(Y)^2 <= C(V):=109^2*V^2*(A^2+2*A*m+1)`.

The original target obeys `|I0|<=Mbar=2*(e+1)^2*(exp(.5)+mbar_delta)`, where the known envelope

`mbar_delta = exp(delta/2 + 1/(2*(1-2*delta)))/sqrt(1-2*delta)`

follows from `z*sqrt(1+z^2)<=z^2+.5` and a Gaussian integral. For nonnegative pi summing to one, `|B0|<=Mbar+|J|+sqrt(C(V))`. Accordingly the conditional MSE of the implemented-form estimator has the bound

`(Mbar+|J|+sqrt(C(V)))^2 + C(V)/N_E`.

A finite inverse residual changes D by its certified log-response-to-index tolerance if that tolerance is available; floating-point roundoff is not formally enclosed here. This bound is uniform in the true original-class profiles and allowed delta, conditional on the finite constructed c,J. The ridge and finite fitted parameter bounds prevent infinite c, but the displayed constant can be enormous and does not tend to zero. It absorbs omitted directions and nonlinear effects very crudely. It supplies no justification for the practical tuning proxy or calibrated intervals and is **not a substantive new rate theorem**.

## 6. Attempt at the requested sharper result — OPEN, gaps explicit

Condition on a learned finite space and write the actual displacement, when meaningful in X, as `h=v+w`, with v in that space. A path expansion formally decomposes B into `R(v)`, the omitted-direction residual `R(w)`, and a nonlinear remainder. Numerical b,G and score/target evaluation add another error. For any justified bounds on these terms, Cauchy–Schwarz converts their squares to squared-target risk units. What is missing is justification of the bounds from this finite source experiment.

| Requested term | Exact object / units | Current status |
|---|---|---|
| Nuisance radius r_N | X norm of the relevant projected displacement; squared radius times squared residual dual norm | No shrinking data-supported radius proved. The fixed `.02` in selection is heuristic, not r_N |
| Operator error | Squared effect of estimated score operator/target derivative and true-versus-fitted variance transfer | Fitted finite geometry is computed; source-distribution-to-operator control with useful constants is unproved |
| Omitted directions | `|R(w)|^2` plus any omitted path contribution | Search supplies lower witnesses only. No finite-bank approximation bound for the original analytic class |
| Nonlinear error | Square of the path-integrated second derivative, or a direct nonlinear B bound | The exact identity and finite paths are checked; no shrinking uniform bound over a source-compatible set |
| Selection error | Excess actual conditional risk of the selected lambda/space over the best frozen choice | No oracle inequality for the source-score proxy; would need source-compatible risk control or independent valid risk estimates |
| Numerical/tail error | Squared target integration, centering, inverse and G/M perturbation effects; also variance effects in squared-target units | Refinement differences and inverse checks saved. M discrepancy is material; no certified total numerical enclosure |

For example, if a true finite-space radius r and bounds o,l,q on omitted, nonlinear and numerical bias were given, `B^2 <=4*(r^2||R||^2+o^2+l^2+q^2)` is immediate. This **conditional implication is not the proposed learned theorem**, because r,o,l,q are not supplied by the observations. Simply naming them as error terms would hide the central problem.

The original envelopes imply a known but very large global X-diameter: coefficient squared distance at most 1, and each profile's squared weighted norm difference at most `pi*(4^2+4^2+20^2+20^2)`. This gives a nonshrinking radius, not the desired r_N. Historical P7–P8 provide a very conservative source-law-based constructive route for full-class identification/consistency, but do not give the neural sieve approximation or learned-score stability needed here. The old matching known-profile lower/upper rate is a baseline/subclass statement; this pilot does not improve it.

## 7. Source compatibility and target-only estimation

For each independent source stratum, DKW gives probability at least `1-alpha` that all five true CDFs satisfy `sup_y |F_a(y)-Fhat_emp,a(y)|<=sqrt(log(10/alpha)/(2*n_a))`, by a union bound. Continuity permits computing the exact empirical/model sup distance at sorted observations. The full original-class set therefore contains truth with this probability. Optimizing I over that full set would target the contrast directly, without requiring a chosen full latent reconstruction.

The implemented endpoint search uses the same source observations and envelope information but restricts to a finite basis family and approximate optimization. Its feasible points are inner witnesses. The sieve may exclude truth; approximate extrema may miss the full range. Neither the midpoint nor those endpoints form a proved confidence interval. A sharp model-specific nonlinear modulus of the target over the **complete source-law experiment**, plus an attainable optimization/approximation guarantee, remains the most consequential constructive gap.

## 8. Retained complete-experiment lower bound

The original unknown-profile class contains the known-exponential-profile subclass. The historical complete-product Hellinger/testing argument therefore still gives minimax squared target risk at least c*min(1,1/(n*delta^2)) for balanced n per source stratum, uniformly over 0<delta<=.1. All five source laws are included; only the designated informative stratum varies in the lower-bound construction. This is a retained subclass lower bound, not a matching learned unknown-profile upper bound. No Wasserstein-to-testing implication is used. See the preserved unknown_profiles/rate_audit.md for its parameter paths and allocation conditions.
