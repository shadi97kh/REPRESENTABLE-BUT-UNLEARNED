# A fixed-profile coefficient path in the original experiment

> Version 2: proof verification continuation, 2026-09-13. Version 1 already existed when this handoff was resumed. Its correct construction is retained here after rederivation; no stronger theorem or new observation model is claimed. All historical files remain unchanged. The current review and execution record is [commands_and_resources.md](commands_and_resources.md).

**PROVED HERE, with verification in [admissibility_and_score.md](admissibility_and_score.md):** the following profiles and the entire coefficient path belong to the original model. Profiles may depend on the known warp when choosing a worst-case truth, but they are fixed as the testing parameter varies. No profile class is narrowed and no observed action is added.

## Experiment and quantifiers

Fix known \(\theta\in[1/2,2]\) and \(0<\delta\le\delta_*:=1/3{,}200{,}000\). Put
\[
s_\theta(z)=z\sqrt{1+\theta z^2},\qquad
v(z)=\delta s_\theta(z),\qquad h(z)=z+v(z).
\]
At each truth \(\eta=(g_1,g_2,\alpha_1,\beta_1,\alpha_2,\beta_2)\), observe exactly five independent groups \(Y_{ai}=F_a(Z_{ai})\), \(i=1,\ldots,n\), \(a=0,\ldots,4\), \(N=5n\). The \(Z_{ai}\) are independent standard Gaussians and unobserved. The observed response curves are
\[
\begin{aligned}
F_0(z)&=g_1(z)+g_2(hz),\\
F_i(z)&=g_1(z+\alpha_i)+g_2(hz+\beta_i),\quad i=1,2,\\
F_3(z)&=g_1(z+1)+g_2(hz),\\
F_4(z)&=g_1(z)+g_2(hz+1).
\end{aligned}
\]
The warp, Gaussian latent law and private unit anchors are known. There is no paired cross-action observation or observed joint-action outcome.

The original profile class is
\[
\mathcal G=\{g:\mathbb R\to(0,\infty)\text{ real analytic}:
g(0)=1,\ g(-\infty)=0,\
\tfrac12e^x\le g'(x)\le2e^x,\
|g''(x)|,|g'''(x)|\le10e^x\}.
\]
All four coefficients remain unknown in \([.5,1]\). The unchanged target and loss are
\[
I(\eta)=E\Delta_{\alpha_1,\alpha_2}g_1(Z)
+E\Delta_{\beta_1,\beta_2}g_2(hZ),\qquad
\ell(\widehat I,I)=(\widehat I-I)^2,
\]
where \(\Delta_{a,b}g(x)=g(x+a+b)-g(x+a)-g(x+b)+g(x)\).

## Explicit hard profiles

Set \(A=1/10\), \(a_0=3/4\), \(a_2=3/5\), \(b_2=9/10\) and \(t_0=1/20\). Define
\[
L_\delta(v)=\delta\left[
\log\cosh\frac{v+A}{2\delta}
-\log\cosh\frac{v-A}{2\delta}\right],\qquad
u(z)=L_\delta(v(z)),
\]
\[
b_\delta=\left[\int_{-\infty}^0e^{x+u(x-a_0)}\,dx\right]^{-1},
\quad
g_{1,\delta}(x)=b_\delta\int_{-\infty}^x e^{y+u(y-a_0)}\,dy,
\quad g_{2,\delta}(x)=e^x.
\tag{1}
\]
These are mathematical definitions, not numerically evaluated integrals. The symbol \(b_\delta\) is a profile normalizer; \(b_2\) is the second fixed coefficient and is different.

For every allowed warp, these profiles obey stronger global bounds than required:
\[
\tfrac45e^x<g_{1,\delta}'(x)<\tfrac54e^x,\qquad
|g_{1,\delta}''(x)|<2e^x,\qquad
|g_{1,\delta}'''(x)|<3e^x.
\tag{2}
\]
Their proof includes the entire saturation transition, rather than using analyticity as a substitute for derivative control.

## Exact nonlinear path

Freeze \(g_{1,\delta},g_{2,\delta},u,b_\delta\) and the internal profile shift \(x-a_0\) throughout \(|t|\le t_0\). Vary only
\[
\alpha_1(t)=a_0+t,\qquad
D(t)=1-b_\delta(e^t-1),\qquad
\beta_1(t)=a_0+\log D(t),\qquad
\alpha_2(t)=a_2,\quad\beta_2(t)=b_2.
\tag{3}
\]
In particular, the profile definition in (1) is never recentered at \(\alpha_1(t)\).
The exact amplitude identity is
\[
e^{\beta_1(t)}=e^{a_0}[1-b_\delta(e^t-1)],\qquad
\frac{d}{dt}e^{\beta_1(t)}=-b_\delta e^{a_0+t}.
\tag{4}
\]
A merely linear path for \(\beta_1\) would not have this identity.

The bounds proved in the companion document are
\[
\frac{161}{171}\le D(t)\le\frac{19}{18},\qquad
\alpha_1(t)\in[.7,.8],\qquad
\beta_1(t)\in(11/16,13/16).
\tag{5}
\]
Thus the entire path is interior to the original coefficient box, and the logarithm is defined everywhere on it.

Exactly four of the five training laws, \(P_0,P_2,P_3,P_4\), are identical for all \(t\). The sole changing curve is
\[
F_{1,t}(z)=g_{1,\delta}(z+a_0+t)+e^{h(z)+\beta_1(t)}.
\tag{6}
\]
This equality is exact at the population-law level, including all nonlinear terms. It does not assert that four laws are generically uninformative in the full model.

## What is being bounded

At a fixed known warp let \(\Theta\) be the full original truth class and define
\[
\mathcal R_{N,\delta,\theta}
=\inf_{\widehat I}\sup_{\eta\in\Theta}
E_\eta[(\widehat I-I_\eta)^2].
\]
The infimum is over all estimators based on the complete five-group experiment. Restricting the supremum to (1)–(3) proves a lower bound for this original problem. The construction is an adversary used in the proof, not information required by an estimator.

Let \(\mathcal H_\eta=\bigoplus_{a=0}^4L_0^2(P_{\eta,a})\), with norm
\(\|\psi\|_\pi^2=(1/5)\sum_aE_a\psi_a^2\). Retain the saved full admissible score closure and define the information bound
\[
V_{\rm eff}(\eta;\delta,\theta)
=\inf\{\|\psi\|_\pi^2:
\langle\psi,Sq\rangle_\pi=DI_\eta[q]
\text{ for every admissible direction }q\}.
\tag{7}
\]
The saved family theorem supplies at least one finite-norm compatible influence at every truth, so this set is nonempty. This definition does not assert that an efficient estimator has been constructed.

The local result concerns the specific truth \(\eta_{\delta,\theta}(0)\) defined above. The truth varies with \(\delta,\theta\) when taking a worst-case supremum. It is not the old exponential-profile reference and does not bound the variance at every truth.

**OPEN:** exact efficient influences, leading exponential constants, matching finite-sample MSE upper bounds and publication priority. The [retained family theorem](../../warp_family/v1/decision.md), old reference-efficiency question and backend no-go are preserved.
