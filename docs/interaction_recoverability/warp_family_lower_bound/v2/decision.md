# Decision: worst-case exponential efficiency verified — v2

**Mathematical decision: worst-case exponential efficiency established; the finite-experiment lower bound is also complete.**

**PROVED HERE, combined with the preserved family upper bound.** For the original truth class \(\Theta\), every fixed known \(\theta\in[1/2,2]\), and \(0<\delta\le1/3{,}200{,}000\), the explicit admissible truth below satisfies
\[
\boxed{
V_{\rm eff}(\eta_{\delta,\theta}(0);\delta,\theta)
\ge \frac{e^{11/5}}{32000}\exp\!\left\{\frac1{640\delta}\right\}.
}
\]
The saved uniform constructed-influence upper bound then yields
\[
\boxed{
\log\left(1+\sup_{\eta\in\Theta}
V_{\rm eff}(\eta;\delta,\theta)\right)=\Theta(1/\delta)
}
\]
as \(\delta\downarrow0\), uniformly in \(\theta\). The lower and upper leading exponential constants are unmatched.

A completed [v1](../v1/decision.md) already existed when this handoff was resumed. This fresh version preserves it and records rederivation of the supplied seed, retained proofs and constants, with expanded checks of saturation, density domination, allocation and quantifiers. No failing inequality was found within this bounded static verification. The mathematical conclusion is unchanged from v1. This is self-review, not independent external review or a formal certificate.

## Exact hard profiles and path

**PROVED HERE.** Set \(A=.1,\ a_0=.75,\ a_2=.6,\ b_2=.9,\ t_0=.05\), and
\[
v(z)=\delta z\sqrt{1+\theta z^2},\qquad
L_\delta(v)=\delta\left[
\log\cosh\frac{v+A}{2\delta}-\log\cosh\frac{v-A}{2\delta}\right],
\qquad u=L_\delta\circ v.
\]
Define
\[
b_\delta=\left[\int_{-\infty}^0e^{x+u(x-a_0)}\,dx\right]^{-1},
\quad
g_1(x)=b_\delta\int_{-\infty}^x e^{y+u(y-a_0)}\,dy,\quad
g_2(x)=e^x.
\]
Both profiles belong to the full original analytic class. Global bounds include
\[
\tfrac45e^x<g_1'(x)<\tfrac54e^x,\qquad
|g_1''(x)|<2e^x,\qquad |g_1'''(x)|<3e^x.
\]
Freeze these profiles, their normalizer and internal shift while varying
\[
\alpha_1(t)=a_0+t,\qquad
\beta_1(t)=a_0+\log[1-b_\delta(e^t-1)],\qquad
\alpha_2=a_2,\quad\beta_2=b_2,\quad |t|\le t_0.
\]
The logarithm is well-defined and all coefficients stay in \([.5,1]\). The [construction](construction.md) retains exactly five independent groups, \(n\) per group, \(N=5n\), the known warp and the original interaction target.

## Complete-experiment control and target separation

**PROVED HERE.** Only source 1 changes. The fixed-outcome score is
\[
S_t(z)=z\,w_t(z)-w_t'(z),\qquad
w_t=\frac{\dot F_{1,t}}{F_{1,t}'}
=\frac{\lambda_t(1-e^{d_t})}{\lambda_t+h'e^{d_t}},
\]
where \(\lambda_t=b_\delta e^{a_0+t-\beta_1(t)}\) and
\(d_t(z)=v(z)-u(z+t)\). For \(\epsilon_\delta=e^{-1/(1280\delta)}\),
\[
\|S_t\|_2\le40(\epsilon_\delta+\delta|t|),\qquad
H(P_{1,t},P_{1,0})\le20\epsilon_\delta|t|+10\delta t^2,
\]
\[
I_\delta'(t)\le-k,\qquad k=e^{11/10}/10.
\]
The [admissibility and score proof](admissibility_and_score.md) verifies positive common support, the density Jacobian and integrable square-root-density derivatives along the whole path. The target margin uses analytic inequalities without quadrature.

**STANDARD RESULT APPLIED.** The complete weighted score norm is \(\|S\|_\pi^2=\|S_0\|_2^2/5\). Cauchy–Schwarz against every full-model compatible influence gives the local variance bound above.

## Separate finite-sample conclusion

**PROVED HERE / STANDARD RESULT APPLIED.** For every integer \(n\ge1\), \(N=5n\), and every allowed \(\delta,\theta\),
\[
\boxed{
\inf_{\widehat I}\sup_{\eta\in\Theta}
E_\eta(\widehat I-I_\eta)^2
\ge \frac{e^{11/5}}{2{,}048{,}000}
\min\left\{1,\frac{e^{1/(640\delta)}}N,
\frac1{\delta\sqrt N}\right\}.
}
\]
The [lower-bound proof](lower_bound.md) gives the stronger constant-specific form. Its choice
\[
t_N=\min\left\{\frac1{20},
\frac1{80\sqrt n\,\epsilon_\delta},
\frac1{\sqrt{40\delta\sqrt n}}\right\}
\]
controls both Hellinger terms. The other four groups have affinity one; the full experiment has total variation at most \(1/2\). Target separation and ordinary squared-loss testing finish the proof with no asymptotic remainder.

The \(\delta t^2\) term produces an intermediate lower-bound scale \(1/(\delta\sqrt N)\). Neither its sharpness nor a true minimax transition is established. An exponential local variance constant does not, by itself, prove exponential sample requirements for fixed accuracy.

## Remaining scope and actual work

**OPEN.** Matching finite-sample MSE upper bounds, exact efficient influences, estimator efficiency, leading exponential constants and publication significance remain unresolved. The hard truth varies with the known warp. Nothing here determines efficiency at the historical \(\delta=.1\) exponential reference or lower-bounds every truth.

The [prior comparison](consequences_and_prior_overlap.md) attributes the information framework to Trabs and uses Heinrich–Kahn only as an adjacent comparison of pointwise and uniform recovery. The possible contribution remains an explicit original-class obstruction in this five-law model together with retained attainment.

All six requested documents are present in this fresh version: [construction](construction.md), [admissibility and score](admissibility_and_score.md), [lower bound](lower_bound.md), [consequences and prior overlap](consequences_and_prior_overlap.md), this decision, and [commands and resources](commands_and_resources.md). Actual work comprised file inspection, analytic rederivation, primary-source retrieval and document authoring. No numerical jobs, samples, fitting, optimization, Gram calculations, benchmarks, GPU/backend work, allowance wrapper, ledger append or sub-agent work occurred. Inspection and authoring overhead was nonzero and was not comprehensively CPU-metered.

The historical family theorem, lower-bound v1, external-review v2, original efficiency audit, failed pilot, predictions/configurations, ledgers and backend no-go remain intact.
