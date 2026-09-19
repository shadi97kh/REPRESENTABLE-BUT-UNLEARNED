# Decision: family attainment established

**Mathematical decision: family attainment established, on the explicit sufficient small-separation interval below. Publication significance remains open.**

This is a proof-only continuation with static internal cross-checks. No numerical evidence was generated. The new documents preserve every historical theorem, negative result, review packet, ledger and implementation decision.

## Strongest completed theorem and all quantifiers

**PROVED HERE.** For every fixed known
\[
\theta\in[1/2,2],\qquad 0<\delta\le1/3{,}200{,}000,
\qquad h(x)=x+\delta x\sqrt{1+\theta x^2},
\]
and every fixed pair of real-analytic profiles satisfying
\[
g_j(0)=1,\quad g_j(-\infty)=0,\quad
\tfrac12e^x\le g_j'(x)\le2e^x,\quad
|g_j''(x)|,|g_j'''(x)|\le10e^x,
\]
with every fixed unknown \((\alpha_1,\beta_1,\alpha_2,\beta_2)\in[.5,1]^4\), the five independent source groups of \(n\) observations each support the explicit estimator in [attainment.md](attainment.md). With \(N=5n\),
\[
\widehat I_N-I_\eta
=\frac1N\sum_{a=0}^4\sum_{i=1}^n\psi_{\eta,a}(Y_{ai})
+o_{P_\eta}(N^{-1/2}),\qquad
E_a\psi_{\eta,a}=0,\quad E_a\psi_{\eta,a}^2<\infty.
\]
The resulting pointwise normal limit has variance
\(V_\eta=\sum_a(1/5)E_a\psi_{\eta,a}^2\).
The influence includes all five source components, both unknown coefficient pairs and combined reference-source covariance. It represents the derivative on the full closure of actual admissible weighted-\(C^3\) profile/coefficient path scores.

Fixed-path regularity holds along each such two-sided DQM \(t/\sqrt N\) path; constrained boundary truths admit only paths staying in the original model. An unconstrained formal direction at a boundary is not silently added to that class.

The estimator uses observed polygonal quantiles, finite commutator reconstruction, fixed public transport charts, public projected calibration and signed target integration. True profiles, their derivatives, true coefficients, Jacobians and influence functions appear only in the proof. A certified finite arithmetic realization is specified for refinable observed inputs and known primitives. No executable backend or fixed-resolution-data guarantee is supplied.

## What extends beyond the calibrated example

**PROVED HERE.** The verified commutator is
\[
C=T^{-1}\circ(\,T(\,\cdot+1\,)-1),\qquad
E(F;x)=g_1(Cx)-g_1(x).
\]
On \(J=[-2,0]\), the general warp lemma gives
\[
C_\delta(x)=x-\delta D_s(x)+r_\delta(x),\quad
\|r_\delta\|_\infty\le12M^2\delta^2,\quad
\|r_\delta'\|_\infty\le246M^2\delta^2,
\]
with invariant \(J\) and
\(0<C'\le1-\kappa\delta/2\) when
\(\delta\le\min\{1/M,\kappa/(500M^2)\}\).
For this warp family, \(M=10,\kappa=1/64\) are verified sufficient choices. All six central source queries are in \([-2,2]\). This proof does **not** cover \(\delta=.1\) and does not claim the sufficient restriction is necessary.

Unlike the old near-reference calibration, the new known radius
\[
R^2=\frac{\sqrt{1+36\theta/\delta^2}-1}{2\theta},
\quad h(R)-R=3,
\]
gives a globally separated scaled population coefficient map over the enlarged search box:
\[
\|B(c)-B(d)\|_\infty\ge\gamma\|c-d\|_\infty,\qquad
\gamma=\tfrac12e^{.4}-2e^{-1.9}>0.
\]
This proof holds for the **full original profile envelope**, without a profile radius shrinking with \(\delta\). Enlarging the computational coefficient box to \([.4,1.1]^2\) does not alter the original truth class.

These arguments, the explicit simultaneous empirical remainder bounds and the new influence constitute the model-specific extension. The separate general \(C^3\) commutator lemma does not by itself prove a statistical theorem for every warp satisfying only its local hypotheses: global moments, calibration and target tails were verified for the stated square-root family.

## Quantitative conditioning and its limits

**PROVED HERE.** [Conditioning](conditioning.md), (4)–(10), gives the fully explicit elementary bound
\[
\|\psi_\eta\|_\pi\le
\Psi(\delta,\theta)=D_t+8D_c+B_I Z_R/\gamma,
\quad
D_c=\frac{6\sqrt2 C_\kappa}{1-\sqrt{1-\kappa\delta/2}},
\]
\[
Z_R=e^R\{2[2Q(2)+D_c]+(4\lceil R+4\rceil+2)Q(R+3)\}.
\]
The linked document defines every remaining public constant and separates central inversion, transports, tail quantiles, row scaling and target integration. In particular
\[
\log(1+\Psi)\le
\frac{3}{4\delta\sqrt\theta}
+O(\delta^{-1/2})+O(\log(1/\delta)).
\]
Thus the sufficient constructed-variance bound grows at most as
\[
V_\eta\le\Psi^2
\le\exp\{3/(2\delta\sqrt\theta)+
O(\delta^{-1/2})+O(\log(1/\delta))\}.
\]
The constants in these asymptotic upper orders can be chosen uniformly in \(\theta\). The explicit \(\Psi\) is available without \(O\)-notation.

The central geometric loss is polynomial. The extreme calibration queries supply the exponential sufficient loss. This is not an exponential impossibility result or an estimate of efficient variance. No practical sample-size regime is demonstrated.

**PROVED HERE / STANDARD RESULT APPLIED.** A rederived exponential-profile coefficient subfamily yields a complete-five-law squared-risk lower bound of order
\[
\inf_{\widehat I}\sup_\eta E_\eta(\widehat I-I_\eta)^2
\gtrsim \min\{1,(N\delta^2)^{-1}\},
\]
with explicit constants in conditioning (15), uniformly in the stated \(\theta\) range. The full nonlinear likelihood path is controlled; the other four source laws are identical along this particular alternative. This is a global subclass lower bound. It is not a pointwise lower bound at every truth.

The polynomial risk lower order and exponential constructed-CLT variance upper order are unmatched. The latter is not an MSE upper bound. The original ordinary squared-error objective is retained, but sharp squared-risk attainment remains open.

## Closest equivalences and significance

**STANDARD RESULT APPLIED.** [Novelty](novelty.md) compares the actual requested primary results. The central series is a normalized cohomological/Neumann inverse; quantile inversion, L-statistics, minimum-distance calibration, preconditioning, adjoint-score regularity and functional recovery without full nuisance recovery are established tools.

Trabs supplies the information/adjoint framework, Kaji the relevant inverse-quantile and integrated-statistic machinery, and Yao–de la Llave the cohomological-series comparison. Bennett and collaborators establish the broader functional-without-full-nuisance principle. DExtrI's main and broader appendix results require their actual observation and rigidity assumptions; five isolated known-warp source laws are not their coordinate-interval experiment.

The concrete addition is a family-level, full-original-profile, four-unknown-coefficient attaining construction in this particular experiment, with its separation costs exposed. Equally informed direct-target, inverse-functional, likelihood and confidence-set controls may use every identity and calibration improvement. No comparative advantage over them is proved.

**OPEN — consequential remaining objection:** is the extreme-calibration exponential loss avoidable, so that the original class admits a sharp or useful separation-dependent attaining result? The present proof does not answer this. In particular the scaled quantile remainder condition displayed in conditioning Section 5 does not follow merely from \(N\delta_N^2\to\infty\).

A valid pointwise family extension addresses the narrow-calibration objection, but does not establish publication priority, a substantive new general ML contribution, efficient recovery, or biological siRNA validity. The historical backend no-go and unresolved-efficiency decisions remain in force. This conclusion triggers no additional backend, training, finite-dictionary or packaging cycle.

## Completed work

The seven new documents are [problem](problem.md), [commutator and calibration](commutator_and_calibration.md), [attainment](attainment.md), [conditioning](conditioning.md), [prior comparison](novelty.md), this decision, and [commands and resources](commands_and_resources.md).

Static internal reviews checked the commutator constants, full-class calibration, tail transfer, actual influence, explicit conditioning and arithmetic specification. They found no remaining concrete contradiction in these bounded checks. This is self-review with separate agents in the same system, not independent external mathematical review. The computation-free analytic conclusions must be evaluated on their proofs.
