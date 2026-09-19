# Decision: worst-case exponential efficiency established

**PROVED HERE.** For every fixed known \(\theta\in[1/2,2]\) and \(0<\delta\le1/3{,}200{,}000\), an explicit truth in the original analytic profile and coefficient class satisfies
\[
\boxed{
V_{\rm eff}(\eta_{\delta,\theta}(0);\delta,\theta)
\ge\frac{e^{11/5}}{32000}
\exp\!\left\{\frac1{640\delta}\right\}.}
\tag{1}
\]
Combined with the preserved uniform constructed-influence upper bound, this proves
\[
\boxed{
\log\!\left(1+\sup_{\eta\in\Theta}
V_{\rm eff}(\eta;\delta,\theta)\right)=\Theta(1/\delta)
}
\tag{2}
\]
as \(\delta\downarrow0\), uniformly for \(\theta\in[1/2,2]\).
The supremum is over the full original profile/coefficient class at the fixed known warp. This is sharp logarithmic order, not a matched leading exponential constant.

## Exact construction and validity

The [construction](construction.md) uses \(A=.1,a_0=.75,a_2=.6,b_2=.9,t_0=.05\):
\[
L_\delta(v)=\delta\left[
\log\cosh\frac{v+A}{2\delta}-\log\cosh\frac{v-A}{2\delta}\right],
\quad u_\delta(z)=L_\delta(\delta s_\theta(z)),
\]
\[
b_\delta=\left[\int_{-\infty}^0e^{x+u_\delta(x-a_0)}dx\right]^{-1},
\quad
g_{1,\delta}(x)=b_\delta\int_{-\infty}^xe^{y+u_\delta(y-a_0)}dy,
\quad g_{2,\delta}(x)=e^x.
\]
The [admissibility proof](admissibility_and_score.md) verifies real analyticity on the entire real line, positivity, normalization, the negative-infinity boundary and all original derivative envelopes globally. In fact
\[
\tfrac45e^x<g_1'<\tfrac54e^x,\qquad |g_1''|<2e^x,\qquad |g_1'''|<3e^x.
\]
These bounds include the saturation transition. No \(\delta\)-shrinking profile neighborhood is imposed.

Keep these profiles and their internal shift fixed while varying
\[
\alpha_1(t)=a_0+t,\quad
\beta_1(t)=a_0+\log[1-b_\delta(e^t-1)],\quad
\alpha_2=a_2,\quad\beta_2=b_2,\qquad |t|\le t_0.
\]
The entire path lies inside the original coefficient box. Sources \(0,2,3,4\) are exactly unchanged. Only source 1 moves.

## Complete likelihood and target control

**PROVED HERE.** The actual fixed-outcome density score retains its Jacobian:
\[
S_t(z)=z\,w_t(z)-w_t'(z),\qquad
w_t=\dot F_{1,t}/F_{1,t}'.
\]
With \(\epsilon_\delta=e^{-1/(1280\delta)}\), the proved uniform bounds are
\[
\|S_t\|_{L^2(\phi)}\le40(\epsilon_\delta+\delta|t|),
\]
\[
H(P_{1,t},P_{1,0})
\le20\epsilon_\delta|t|+10\delta t^2,\qquad
I_\delta'(t)\le-\frac{e^{11/10}}{10}.
\tag{3}
\]
Positive common support and a common-coordinate dominated square-root-density derivative justify integration along the whole nonlinear path. No unproved nonlinear remainder or Wasserstein-to-testing inference is used.

The score norm at \(t=0\) is strictly positive and exponentially small. Its full-experiment norm has the allocation factor \(1/5\). The standard information inequality gives (1) for every full-model compatible influence, without a finite tangent dictionary.

## Separate finite-sample conclusion

**PROVED HERE / STANDARD RESULT APPLIED.** For every integer \(n\ge1\), \(N=5n\), throughout the same \(\delta,\theta\) range,
\[
\boxed{
\inf_{\widehat I}\sup_{\eta\in\Theta}
E_\eta[(\widehat I-I_\eta)^2]
\ge\frac{e^{11/5}}{2{,}048{,}000}
\min\left\{
1,\frac{e^{1/(640\delta)}}N,\frac1{\delta\sqrt N}
\right\}.}
\tag{4}
\]
The stronger constant-specific form, exact test separation, five-group product affinity and squared-loss reduction are in [lower_bound.md](lower_bound.md). Both terms in the Hellinger bound determine the separation; the nonlinear term has not been discarded.

This lower bound does not imply exponential sample requirements for fixed accuracy. The intermediate scale \(1/(\delta\sqrt N)\) can govern before the much later local-information regime, with crossover of order \(\delta^2e^{1/(320\delta)}\) for this bound. No matching finite-sample upper bound is proved.

## What changes and what remains open

The previous conjectural possibility of a polynomial uniform bound on **worst-case regular local variance over the original class** is ruled out. Exponential worst-case local conditioning is not solely an artifact of the earlier estimator. This is a substantive additional mathematical obstruction within the investigated model.

The hard truth depends on \(\delta,\theta\). The result does not lower-bound efficient variance at every truth, does not determine the historical exponential reference at \(\delta=.1\), and does not prove the retained estimator efficient. The earlier family-attainment theorem remains intact: finite pointwise regular variance at each fixed allowed warp is consistent with an exponentially diverging worst-case bound.

**OPEN — principal remaining gap:** matching finite-sample minimax MSE upper bounds, including the nonlinear regime and separation dependence. Exact efficient influences, leading exponential constants, practical error guarantees, publication priority and biological validity also remain unresolved.

[Prior overlap](consequences_and_prior_overlap.md) identifies the standard submodel information reduction (Trabs) and the established distinction between pointwise and uniform rates (Heinrich–Kahn). Their finite-mixture rates do not apply automatically here. Smooth saturation, higher-order cancellation and Hellinger testing are established tools. The concrete contribution candidate is this explicit original-class five-law obstruction together with retained attainment; no new general ML principle or publication priority is claimed.

The historical reference-efficiency question, failed learned-target pilot, frozen predictions/configurations, external-review v2 and backend no-go remain unchanged. No estimator, backend, training, dictionary or packaging cycle is initiated.

## Review and execution status

Completed files: [construction](construction.md), [admissibility and score](admissibility_and_score.md), [lower bounds](lower_bound.md), [consequences and prior overlap](consequences_and_prior_overlap.md), this decision, and [commands and resources](commands_and_resources.md).

The derivations were checked through static self-review with separate agents inside the same system. This is not independent external mathematical review or a formal proof-assistant certificate. No numerical job, sample, fit, GPU operation or ledger write was performed.
