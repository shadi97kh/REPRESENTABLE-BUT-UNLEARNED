# Local information and finite-experiment squared-risk lower bounds

The parameter ranges, original model and hard path are fixed in [construction.md](construction.md). Its full global admissibility, actual density score and nonlinear path regularity are proved in [admissibility_and_score.md](admissibility_and_score.md). This document uses no numerical output.

Put
\[
c_0=1/1280,\quad
\epsilon_\delta=e^{-c_0/\delta},\quad
k=e^{11/10}/10,\quad N=5n,\quad n\ge1.
\]
The verified inputs, uniformly throughout \(|t|\le t_0=1/20\), are
\[
\|S_t\|_{L^2(\phi)}\le40(\epsilon_\delta+\delta|t|),
\quad
H(P_{1,t},P_{1,0})\le20\epsilon_\delta|t|+10\delta t^2,
\]
\[
I_\delta'(t)\le-k<0.
\tag{1}
\]
The other four source laws are exactly constant along this path.

## 1. Exponential efficient-variance obstruction at the hard truth

**STANDARD RESULT APPLIED / PROVED HERE for the model-specific inputs.**
At the truth \(\eta_{\delta,\theta}(0)\), the coefficient direction is
\[
(d_{\alpha_1},d_{\beta_1},d_{\alpha_2},d_{\beta_2})
=(1,-b_\delta,0,0).
\]
Its profiles do not vary. It is a genuine two-sided admissible direction, since the entire path is interior to the coefficient box and the profiles have strict derivative margins.

The full five-source score vector has only source 1 nonzero, and its norm is
\[
\|S\|_\pi^2=\tfrac15 E_\phi S_0(Z)^2
\le320\epsilon_\delta^2.
\tag{2}
\]
The subscript 0 on \(S_0\) means the path parameter \(t=0\), not source law zero. Source law zero has zero score in this direction.

Every full-model compatible influence satisfies
\(\langle\psi,S\rangle_\pi=I_\delta'(0)\).
Cauchy–Schwarz therefore gives
\[
\|\psi\|_\pi^2\ge
\frac{|I_\delta'(0)|^2}{\|S\|_\pi^2}
\ge\frac{k^2}{320\epsilon_\delta^2}.
\]
Taking the infimum over compatible influences yields
\[
\boxed{
V_{\rm eff}(\eta_{\delta,\theta}(0);\delta,\theta)
\ge \frac{e^{11/5}}{32000}
\exp\!\left\{\frac{1}{640\delta}\right\}.}
\tag{3}
\]
This holds for every allowed \(\theta,\delta\). It is a bound at the explicit delta-dependent truth, not at every original-class truth and not at the historical exponential reference.

The nonempty compatible-influence set follows from the retained family theorem, whose assumptions include this verified profile pair. Thus (3) is a finite but exponentially diverging information obstruction, compatible with existence of a pointwise regular estimator. The score norm is strictly positive at each positive \(\delta\), as proved in the companion document. No zero-information assertion or finite-dictionary argument is used.

The information inequality and its interpretation through the full score closure are standard; the concrete exponentially small score with a nonvanishing target derivative is the new model-specific construction.

## 2. All five training groups and nonlinear Hellinger control

**PROVED HERE / STANDARD RESULT APPLIED.** Denote the complete training distribution by
\[
\mathbb P_t=\bigotimes_{a=0}^4P_{a,t}^{\otimes n}.
\]
For the specified path \(P_{a,t}=P_{a,0}\) exactly when \(a\in\{0,2,3,4\}\). Thus their affinities equal one, and only the \(n\) observations from source 1 contribute to the product affinity.

With \(H^2(P,Q)=\int(\sqrt p-\sqrt q)^2\), affinity is
\(\rho(P,Q)=1-H(P,Q)^2/2\). Hence
\[
H(\mathbb P_t,\mathbb P_0)^2
=2[1-\rho(P_{1,t},P_{1,0})^n]
\le nH(P_{1,t},P_{1,0})^2.
\]
Also
\(\operatorname{TV}(P,Q)\le H(P,Q)\), by Cauchy–Schwarz applied to
\((\sqrt p-\sqrt q)(\sqrt p+\sqrt q)\) and the factor \(1/2\) in TV.
Together with the actual nonlinear segment bound (1),
\[
\operatorname{TV}(\mathbb P_t,\mathbb P_0)
\le \sqrt n\,[20\epsilon_\delta|t|+10\delta t^2].
\tag{4}
\]
This is a finite-experiment result. It does not replace Hellinger distance with Wasserstein proximity, nor use a pointwise Taylor remainder without controlling the segment.

Choose the positive path separation
\[
t_N=\min\left\{
\frac1{20},\
\frac1{80\sqrt n\,\epsilon_\delta},\
\frac1{\sqrt{40\delta\sqrt n}}
\right\}.
\tag{5}
\]
The first term keeps the whole segment admissible. The second makes the first contribution in (4) at most \(1/4\); the third makes the second contribution at most \(1/4\). Therefore
\(\operatorname{TV}(\mathbb P_{t_N},\mathbb P_0)\le1/2\).
The two target values differ by at least \(kt_N\).

## 3. Explicit finite-sample minimax lower bound

**STANDARD RESULT APPLIED.** For any estimator choose the test that selects the closer of the two target values. A wrong selection entails squared error at least one quarter of their squared separation. The sum of the two test errors is at least \(1-\operatorname{TV}\), and the maximum of two risks is at least half their sum. Consequently
\[
\max_{t\in\{0,t_N\}}E_t(\widehat I-I_\delta(t))^2
\ge\frac{k^2t_N^2}{8}
[1-\operatorname{TV}(\mathbb P_{t_N},\mathbb P_0)]
\ge\frac{k^2t_N^2}{16}.
\]
The supremum over the original class is at least this two-point maximum, for every estimator using all five observed groups. In total-\(N\) notation, \(n=N/5\), so
\[
\boxed{
\mathcal R_{N,\delta,\theta}\ge
\frac{k^2}{16}
\min\left\{
\frac1{400},\
\frac{1}{1280N\epsilon_\delta^2},\
\frac{\sqrt5}{40\delta\sqrt N}
\right\}.}
\tag{6}
\]
This holds for every integer \(n\ge1\), every \(\theta\in[1/2,2]\), and every \(0<\delta\le1/3{,}200{,}000\), with no asymptotic remainder.

Since \(1/1280\) is no larger than the other two scalar constants inside the minimum, a simpler, weaker bound is
\[
\boxed{
\mathcal R_{N,\delta,\theta}\ge
\frac{e^{11/5}}{2{,}048{,}000}
\min\left\{
1,\
\frac{\exp\{1/(640\delta)\}}{N},\
\frac1{\delta\sqrt N}
\right\}.}
\tag{7}
\]
Every constant is explicit and independent of \(\theta,\delta,N\) over the declared ranges. The hard profiles may vary with known \(\delta,\theta\) in a worst-case argument, but the same pair is held fixed at both testing endpoints and along their connecting path.

No oracle data are used to make two laws close. Indeed, even disclosing the fixed hard profiles to an estimator would not change the two endpoint distributions or the testing inequality for this subexperiment. The result for the original unknown-profile problem requires no such disclosure.

## 4. Relation to the saved upper bound

**PROVED HERE, using the preserved family-attainment theorem.** At each original-class truth the saved construction has compatible influence \(\psi_\eta\) with
\[
V_{\rm eff}(\eta;\delta,\theta)\le\|\psi_\eta\|_\pi^2
\le\Psi(\delta,\theta)^2,
\]
where \(\Psi\) is completely explicit in
[warp_family/v1/conditioning.md](../../warp_family/v1/conditioning.md), (4)–(10), and uniform over the original profile and coefficient class at each known warp. Therefore
\[
\frac{e^{11/5}}{32000}e^{1/(640\delta)}
\le \sup_{\eta\in\Theta}V_{\rm eff}(\eta;\delta,\theta)
\le\Psi(\delta,\theta)^2.
\tag{8}
\]
The saved bound further gives
\[
\log(1+\Psi^2)\le
\frac{3}{2\delta\sqrt\theta}
+O(\delta^{-1/2})+O(\log(1/\delta)),
\]
with constants uniform in \(\theta\in[1/2,2]\). The positive lower exponential in (8), after absorbing its fixed prefactor, gives
\[
\boxed{
\log\!\left(1+\sup_{\eta\in\Theta}
V_{\rm eff}(\eta;\delta,\theta)\right)
=\Theta(1/\delta)\quad(\delta\downarrow0),}
\tag{9}
\]
uniformly in the stated \(\theta\) range. More precisely there exist constants \(0<c<C<\infty\) and \(\delta_0>0\), independent of \(\theta\), giving the corresponding two inequalities for every \(0<\delta\le\min(\delta_0,\delta_*)\).

This is a matched **logarithmic order of worst-case efficient variance**, not matching leading exponential constants. It does not assert that the constructed influence is efficient, that its pointwise CLT is uniform over the class, or that it gives a finite-sample MSE upper bound matching (6) or (7).

## 5. The nonlinear intermediate regime

The term \(\delta t^2\) must be retained in the proved Hellinger upper bound (4); its sharpness as an intrinsic Hellinger order has not been established. The sufficient separation selected in (5) is, up to constants, the minimum of these three scales:
\[
1,\quad (\sqrt N\,\epsilon_\delta)^{-1},\quad
N^{-1/4}\delta^{-1/2}.
\]
Their squared scales produce the minimum in (7). The two nonconstant risk terms cross at
\[
N\asymp\delta^2\epsilon_\delta^{-4}
=\delta^2 e^{1/(320\delta)}.
\tag{10}
\]
For intermediate counts, schematically between a constant multiple of \(\delta^{-2}\) and this crossover, the nonlinear term can govern the lower bound:
squared risk of order at least \(1/(\delta\sqrt N)\), or root-mean-squared error of order at least \(\delta^{-1/2}N^{-1/4}\).
These are regimes of the proved lower bound; no matching estimator has been supplied.

At each fixed \(\delta>0\), eventually the local information term proportional to \(\epsilon_\delta^{-2}/N\) governs this bound. The onset can depend exponentially on inverse separation. Consequently the exponential information constant alone does not establish exponential sample requirements for any fixed error tolerance. The intermediate lower-bound term already declines at polynomial counts in \(1/\delta\); a lower bound alone supplies no sufficient sample complexity.

**OPEN:** matching finite-sample MSE upper bounds, leading exponential constants, efficiency at the old fixed reference, efficient estimation and useful practical performance. The precise worst-case efficiency question in (9) is resolved within the stated model; these other questions are not.
