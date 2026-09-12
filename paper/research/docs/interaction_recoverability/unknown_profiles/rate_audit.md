> Publication copy of `docs/interaction_recoverability/unknown_profiles/rate_audit.md`. Mathematical and result content is retained; execution material is omitted. This dated record retains its original conclusions and limitations, including conclusions superseded by later stages. Original SHA-256: `e7437820236bea6892388b10fc76e2932fc1d906559e907049eb383065186050`.

# Independent audit of the historical rate

The historical documents and decision are preserved. This audit supports their rate order; it makes allocation and boundary qualifications explicit. It does not transfer that rate to unknown profiles.

| Claim | Audit outcome and mathematical status |
|---|---|
| Balanced known-profile rate | Supported — **STANDARD RESULT APPLIED** after checking the density, path and quantile arguments below. |
| Meaning of n | Clarified — **PROVED HERE**: n per stratum, total N=5n. Total N alone is insufficient for arbitrary allocations. |
| Uniformity as delta tends to zero | Supported for the separate known-delta classes, 0<delta<=.1; delta=0 excluded. **PROVED HERE**. |
| Fixed T/Q lower bound implies the high-information rate | Would be false in isolation. The actual proof additionally supplies a shrinking amplitude path, which suffices. **STANDARD RESULT APPLIED**. |
| Protected contrast has O(1/n) MSE | Supported only in the exponential-profile class; it is an observed-mean identity. **PROVED HERE**. |
| Unknown-profile rate or optimal constants | Not established by the old result. **OPEN**. |
| Numerical stability proves statistical indistinguishability | Not claimed or used. The likelihood/Hellinger proof concerns the complete experiment. **STANDARD RESULT APPLIED**. |

**PROVED HERE — experiment and range.** The old class is
\(F_a=e^{z+a_3}u_1^{a_1}u_2^{a_2}+e^{h_\delta(z)+a_4}v_1^{a_1}v_2^{a_2}\), with four unknown coefficients in \([m,e]\), \(m=\sqrt e\). Gaussian noise law, profiles exp, warps, delta, supports, unit private loadings and zero additive terms are known. All observations are independent, with unobserved within-response common latent noise. The target is
\(I=m_0(u_1-1)(u_2-1)+m_\delta(v_1-1)(v_2-1)\).
Its exact range is \([(m-1)^2(m_0+m_\delta),(e-1)^2(m_0+m_\delta)]\). The loss is ordinary squared error, without normalization. Counts at 0,e3,e4 carry no information about these four unknown coefficients; this statement changes in the extension.

**PROVED HERE — audited score bounds.** For \(F_t=Ae^z+Be^{z+ts(z)}\), differentiate the density at fixed outcome, not fixed latent coordinate. With \(u=\dot F/F'\), the score pulled back to latent space is \(zu-u'\). Supports are all \((0,\infty)\); there are no missing atoms or moving endpoints. For \(r=(B/A)e^{ts}\), \(u=s/(r^{-1}+1+ts')\). Using \(|s|\le z^2+1/2\), \(s'\le1+2|z|\), \(|s''|\le2\), gives the old polynomial envelope \(P=1.2|z|^3+.3|z|^2+2.6|z|+1.15\) and \(C_0=EP^2=65.6233873491624\).

The amplitude path fixes \((u_1,v_1)=(m,e)\) and varies \((u_2,v_2)=(b,m+e-b)\). To check its derivative bound more explicitly, put \(E=e^{\delta s}\), \(D_b=b+(m+e-b)E(1+\delta s')\). Then
\[
w=(1-E)/D_b,\quad |w|\le\delta|s|/m,
\quad |D_b'/D_b|\le\delta(s'+2),
\quad |E'/D_b|\le\delta s'/m.
\]
Therefore \(|w'|\le(\delta/m)\{s'+.1|s|(s'+2)\}\), as required. Normal score second moments bound the L2 derivative of square-root densities; integrating that derivative establishes the Hellinger bound without assuming numerical quadrature is exact.

**STANDARD RESULT APPLIED — both lower-bound regimes, all strata.** Use \(H^2=\int(\sqrt p-\sqrt q)^2=2-2\rho\), so \(TV\le H\) and affinities multiply. Only stratum e2 changes along the amplitude path. Let \(D=e-m\), \(\kappa=m_0D\), and choose two centered b values separated by
\[
\eta=\min\{D,m/(\delta\sqrt{C_0n_2})\},
\]
with \(\eta=D\) when \(n_2=0\). These parameters stay in the original box. Single-observation \(H\le\eta\delta\sqrt{C_0}/(2m)\), hence complete-experiment \(TV\le1/2\). The target derivative has magnitude at least \(\kappa\). Thresholding any estimator gives risk at least \(\kappa^2\eta^2/16\). At low information the separation stays fixed; at high information it shrinks. Symmetrically vary the coefficients at e1. Consequently, with \(n_\dagger=\min(n_1,n_2)\),
\[
R_{\boldsymbol n,\delta}\ge c\min\{1,(n_\dagger\delta^2)^{-1}\},
\]
interpreting the right side as a positive constant if \(n_\dagger=0\). This is a bound for every procedure using all strata. Neither Wasserstein proximity nor a fixed alternative alone proves it.

**STANDARD RESULT APPLIED — upper bound and uniform constants.** The central quantile matrix determinant is exactly \(2\sinh(\delta s(r))\), \(r=\Phi^{-1}(.75)\); its bounded entries imply inverse norm at most \(C/\delta\). The response quantile derivative \(F'(z)/\phi(z)\) is uniformly bounded for latent \(z\in[\Phi^{-1}(.125),\Phi^{-1}(.875)]\), uniformly over coefficients and delta. Separate DKW inequalities give
\(P(D_n>t)\le2e^{-2n_1t^2}+2e^{-2n_2t^2}\le4e^{-2n_\dagger t^2}\).
On \(D_n\le1/8\), inversion followed by coordinate clipping gives the claimed Lipschitz error. Off this event the parameter box gives a deterministic bounded error. Integrating the DKW tail and absorbing \(e^{-n_\dagger/32}\) yields the matching upper order. If one informative stratum has zero observations, use the target-range midpoint for the constant upper bound.

Thus
\[
R_{\boldsymbol n,\delta}\asymp\min\{1,(n_\dagger\delta^2)^{-1}\}.
\]
Equivalently one may use \(\delta^{-2}(1/n_1+1/n_2)\), within a factor two before truncation. If both informative allocations are at least \(\pi_*N\), constants can depend on \(\pi_*^{-1}\). No bound with an allocation-free N is asserted. All joint sequences of positive delta and sample counts are covered; the inverse penalty is explicit. Second moments are uniform because \(2\delta<1/2\). The fifth moment at the endpoint is infinite, but the proof does not require it. Computing \(m_\delta\) approximately adds deterministic moment error; it is exact in the statistical theorem.

**PROVED HERE — protected contrast audit.** In the old class only,
\(J=(e-1)(\mu_1-\mu_0)\), so independent sample means have variance
\((e-1)^2\{\operatorname{Var}(Y_1)/n_1+\operatorname{Var}(Y_0)/n_0\}\).
This is standard cancellation of component directions and needs positive n0,n1. Unknown profiles generally do not satisfy \(g(x+1)=eg(x)\); this identity is not carried into the extension.

**NUMERICALLY CHECKED — reused evidence.** The existing four-delta output and 14-test record are retained, not rerun. At delta .0001 they report axis W1 .0004181120591526363, joint W1 1.887095558438503, interaction gap 1.8867185830831366 and axis Hellinger squared 1.2974652568835723e-9. The common limit is 1.886070833180237. These ordinary floating-point diagnostics are not certified enclosures or a proof. See the original [diagnostics](../../../../../runs/interaction_recoverability/prep-20260912-082154/numerics-v2/diagnostics.json). The original optional diagnostic script remains unaudited; its availability is not needed for this proof audit.
