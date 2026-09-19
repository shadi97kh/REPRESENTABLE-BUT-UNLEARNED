# Learned interaction estimator, version 1

Status: implemented and development-tested; frozen pilot results are recorded separately in `results.md`. This is a new authorized experiment, in a separate namespace. The historical known-profile theorem, unknown-profile P1–P8, negative results and disabled interfaces retain their meaning.

## Experiment and learned functions

The experiment is the original [unknown-profile model](../unknown_profiles/model.md): independent observations from five laws indexed by `0,e1,e2,e3,e4`, known standard Gaussian latent law, known `h_delta(z)=z+delta*z*sqrt(1+z*z)`, and known unit private-anchor loadings. Each response's two components share its latent variable; different responses have independent, unobserved latent draws. The target is the expected `e1+e2 - e1 - e2 + 0` contrast on the original positive outcome scale. No joint-action observations enter training.

The four fitted parameters are **index shifts**, ordered in code as `(alpha1,beta1,alpha2,beta2)`. They are the logarithms of the multiplicative coefficients in the older exponential-profile submodel. They are projected onto `[.501,.999]`, a strict subset of the original `[.5,1]` box; this creates a boundary approximation limitation. Both profiles and all four shifts are fitted.

For each profile the width-four neural family is

`g(x)=exp(x)*(1 + .1/4 * sum_j a_j*(tanh(w_j*x+b_j)-tanh(b_j)))`,

with `|a_j|<=1`, `.25<=w_j<=1`, `|b_j|<=3`. Adam steps are followed by exact parameter projection. This additive modulation was chosen in place of the prompt's optional exponential modulation. It is a nonlinear learned profile: amplitudes, frequencies/slopes and locations are all trainable, not fixed oracle profiles.

**PROVED HERE — global admissibility.** Writing `g=exp(x)(1+q)`, tanh derivative bounds give `|q|<=.2`, `|q'|<=.1`, `|q''|<=.2`, `|q'''|<=.2`. Consequently `g'/exp(x)` is in `[.7,1.3]`, `|g''|/exp(x)<=1.6`, and `|g'''|/exp(x)<=2.3`. Also `g(0)=1`, `g(-infinity)=0`, positivity and real analyticity hold. These are global bounds, not conclusions from grid checks. This finite family is substantially smaller than the original analytic class.

The equally informed basis control fits `g=exp(x)*(1+.08*theta^T k(x))`, with each profile's `||theta||_1<=1`. The 16 atoms are six sines (frequencies `.5,1,1.5,2,2.5,3`), seven centered Gaussian bumps (width `.75`, centers `-3,...,3`) and three tanh atoms (width 1, centers `-2,0,2`), each normalized to vanish at zero. The global bounds `|k|<=2`, `|k'|<=3`, `|k''|<=9`, `|k'''|<=27` imply derivative bounds `[.6,1.4]`, `2.36`, `6.2` respectively. Projection onto the L1 ball enforces these restrictions.

## Source objective and numerical inverse

The training objective is the proper population CRPS, `E|Y-F(Z)| - .5 E|F(Z)-F(Z')|`, averaged with the source stratum proportions. Implementation uses 64 deterministic Gauss–Hermite nodes; the monotonicity of F simplifies the pair term to a weighted sum. It is a discretized proper objective, with a saved 64-versus-128 integration diagnostic, not an exact finite-node likelihood. Adam uses learning rate `.03` and 300 steps. Source-selection data choose among checkpoints 100, 200, 300. Pooled controls use the last checkpoint and all observations. There are no target-selected restarts or architecture choices.

Likelihood is implemented for the score/gradient checks, but is not the fitted objective. The monotone inverse uses 52 bisections on `[-24,24]`, with bracket and relative response-residual checks. At fixed outcome y its first derivative is `dz/dtheta = -dotF/F'`. The differentiable inverse attaches precisely that first derivative to the bisection solution. This is sufficient for the first-order optimizers used here; it is not advertised as an exact higher-order implicit differentiator.

## Cross-fitted correction

Each stratum is randomly split into three disjoint parts. In rotation k, fit uses k, selection uses k+1, and correction evaluation uses k+2 (modulo 3). Both the correction space and lambda are chosen before this rotation's evaluation outcomes are accessed. The final estimate averages the three rotation estimates. Their dependence precludes treating their variance estimates as independent.

For a direction v, use the actual fixed-outcome score

`s_v(F(z)) = z*dotF/F' - dotF'/F' + dotF*F''/(F')^2`.

Scores are centered separately under each fitted source law. With `pi_a=n_a,E/N_E`, the fitted Gram G, target derivative b, and direction penalty M give `(G+lambda*M)c=b` on the current space. The estimate is the plug-in **plus** `sum_a pi_a*mean_E psi_a`. The actual stratum variance contribution is `sum_a pi_a^2*Var(psi_a)/n_a,E`; code retains this expression. No correction clipping is applied. Score-null target components remain in the reported residual.

The primary score and target integrals use 96 Legendre nodes on `[-10,10]` against the Gaussian density, renormalized on that interval; 192-node differences are saved. This numerical centering differs from exact full-law centering. Refinement differences are diagnostics, not certified quadrature/tail bounds.

## Independent direction search and controls

The search bank has 36 directions: all four shifts and the 16 analytic atoms multiplied by exp(x) for each profile. Its oscillatory/localized directions extend beyond the fitted network's parameter derivatives. It is **not** claimed to contain the entire neural Jacobian as a subspace, or to exhaust the infinite-dimensional tangent space. The initial fixed correction uses the four shifts plus two atoms for each profile (eight dimensions).

The metric approximates the original weighted Sobolev norm, integrating derivative orders zero through three with weight `exp(-2*x)/(1+x*x)`. After `x=tan(theta)`, 4096 Legendre nodes are used, with a 2048-node comparison. The measured discrepancy is retained. Algebraic statements about unit vectors and spectral penalties refer exactly to this finite numerical M unless explicitly stated otherwise; they are not certified original-X normalizations.

For residual `R=b-Gc`, the adaptive search solves `v proportional to (G+tau*M)^(-1)R`, with `tau=.001`, and M-normalizes it. This solves the finite generalized Rayleigh quotient; it is standard algebra, not a new direction-learning theorem. Accepted nonredundant directions are appended and the correction re-solved, for at most four rounds. The matched random control substitutes an independent Gaussian right-hand side in the same solve and uses the same round limit. The fixed richer control uses all 36 directions. Redundancy tolerance is `1e-9` in M-orthogonalization; a numerically singular whole-bank metric raises an error rather than deleting a target direction silently.

Three lambdas (`1e-5,.001,.1`) are compared using selection-fold scores: a regularized local displacement estimate, the residual's value on that displacement, model variance/N_E, and `.02^2` times the squared whitened residual. This is a **heuristic**, with no justified nuisance radius or oracle selection inequality. It supplies no calibrated confidence interval or full-class risk certificate.

The same neural nuisance fits support the shared plug-in, fixed, adaptive, random and rich controls. The basis family has shared/pooled plug-ins and a rich corrected control. The direct-target basis control searches both target endpoints subject to five DKW source-law constraints on all observations, using 40 optimization steps per endpoint; it returns a midpoint only when both searches find a feasible iterate. Feasibility is checked at every sorted outcome. This is an approximate **inner** search in a restricted sieve, not certified confidence-set inversion over the original class. Failed feasibility and unproved endpoint optimality are reported.

## Admissible finite paths and costs

An enriched direction couples shifts and both profiles. The finite path changes coefficients and adds `t*r_j` to each profile, then evaluates at the changed indices; it includes their nonlinear cross terms. A symmetric path length is bounded by the remaining coefficient-box margin and by `.9*min(.1/5,3.8/65)/max_j ||v_profile,j||_1`. The global derivative margins above prove admissibility. Saved positive and negative finite changes compare exact bias integrals (numerically evaluated) with the linear residual and check selection-data DKW compatibility. A found path is a witness, never a supremum upper bound.

Fully charged per-method costs include the shared nuisance fits and geometry needed by that method. Actual total run cost counts shared work once, includes all endpoint searches, diagnostics and failures, and comes from the process-tree ledger. Summing fully charged comparison costs intentionally overcounts shared work. These costs measure this CPU implementation; three replicates do not establish a general computational scaling result.

## Post-pilot static audit

The candidate nonredundancy precheck uses norm tolerance 1e-7, while the subsequent M-Gram eigensolver uses rank tolerance 1e-9. Passing the former does not guarantee surviving the latter. The 288 stored accepted appends yielded final dimensions 10-12 rather than always 12. These are candidate appends, not 288 certified rank increases; the mismatch is retained as an implementation limitation without post-evaluation refitting. Every corrected arm computes the full 36-direction geometry, so the present cost figures do not demonstrate savings from avoiding that bank.
