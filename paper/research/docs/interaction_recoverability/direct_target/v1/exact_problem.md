> Publication copy of `docs/interaction_recoverability/direct_target/v1/exact_problem.md`. Mathematical and result content is retained; execution material is omitted. This dated record retains its original conclusions and limitations, including conclusions superseded by later stages. Original SHA-256: `e239593b78d32b77ac5883bd184174a18f3b06ea3d344ce14d8ac9d1ae16eadf`.

**Version 1, 2026-09-12. Decision: Unresolved.** This bounded continuation preserves the original class and the completed pilot. Claim labels distinguish proofs from established tools, numerical diagnostics and open steps.

**PROVED HERE — declaration of the experiment and target.** Fix known computable `0 < delta <= .1`, `h(z)=z+delta*z*sqrt(1+z^2)`, and independent unobserved standard Gaussian `Z_ai`. Observe `Y_ai=F_a(Z_ai)` for precisely `a in {0,e1,e2,e3,e4}`, with fixed counts `n_a`, where

`F_a(z)=g1(z+alpha1*a1+alpha2*a2+a3)+g2(h(z)+beta1*a1+beta2*a2+a4)`.

All four shifts lie in `[1/2,1]` and are unknown. Each unknown profile is positive and real analytic, `g(0)=1`, `g(-infinity)=0`, `.5 exp(x)<=g'(x)<=2 exp(x)`, and `|g''(x)|,|g'''(x)|<=10 exp(x)`. The calibrated private shifts are exactly one. There are no observed latent coordinates, cross-action pairings, unknown warps or axis intervals. The two summands of one response share its latent draw; different observations are independent. No joint-action measurements are supplied to an estimator.

The primary target and loss remain

`I12=E[F_(e1+e2)(Z)-F_e1(Z)-F_e2(Z)+F_0(Z)]`, `loss=(Ihat-I12)^2`.

Proof of the basic properties: integrating the derivative envelope from negative infinity gives `.5 exp(x)<=g(x)<=2 exp(x)`. Thus every response is strictly increasing onto `(0,infinity)`, so `q_a(z)=Q_a(Phi(z))=F_a(z)`. The five-law experiment is the product of the five independent stratum experiments. Profile normalization and the coefficient box are unchanged; no analytic closure or neural sieve replaces this class.

Write `m_(p,1)=E exp(pZ)=exp(p^2/2)` and `m_(p,2)=E exp(p h(Z))`. For `2p delta<1`, Gaussian completion of the square gives the known bound

`m_(p,2) <= exp(p delta/2 + p^2/[2(1-2p delta)]) / sqrt(1-2p delta)`.

This follows from `h(z)<=z+delta(z^2+1/2)`. Both first and second moments are uniformly finite for the allowed delta. We use these known envelope integrals, never latent samples. A valid target radius is `Mbar=2(e+1)^2(m_(1,1)+m_(1,2))`, with the displayed upper bound substituted for `m_(1,2)` when a closed-form computable ceiling is desired.

**STANDARD RESULT APPLIED — allocation conventions.** Set `N=sum n_a` and `pi_a=n_a/N`. Population modulus statements assume all `pi_a>0`. The complete Hellinger metric is `d_pi^2=sum_a pi_a H_a^2`, with `H^2(P,Q)=integral(sqrt(p)-sqrt(q))^2`. Balanced `n` means per stratum, so `N=5n`. If `n1` or `n2` is zero, the retained exponential subclass has a positive nonidentification lower risk. If another stratum is absent, the construction here falls back to a full-range answer; no general impossibility theorem for that missing-stratum pattern is inferred.

**PROVED HERE — scope of the attempted construction.** Coefficients are learned from an independent source split using the existing finite-grid calibration argument. The target evaluation stage integrates known signed weights against empirical source quantiles; it does not reconstruct complete profiles. Transient anchor reconstructions at coefficient-calibration queries remain necessary for the available coefficient guarantee. A known-coefficient calculation is a conditional building block, not the final unknown-coefficient result.

Required evidence was inspected: [original model](../../unknown_profiles/model.md), [P1–P8](../../unknown_profiles/proofs.md), [rate audit](../../unknown_profiles/rate_audit.md), learned-target [method](../../learned_target_v1/method.md), [theorem](../../learned_target_v1/theorem.md), [novelty](../../learned_target_v1/novelty.md), [results](../../learned_target_v1/results.md), [review](../../learned_target_v1/review.md), [subsequent adversarial audit](../../learned_target_v1/adversarial_review_20260912.md), source manifests, and actual unknown-profile and learned-target code. The source comparison is in [novelty.md](novelty.md); the sharp unresolved inequality is in [modulus.md](modulus.md).
