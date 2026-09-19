**STANDARD RESULT APPLIED / PROVED HERE — feasible direct finite estimator with a conservative finite-sample bound.** All statements use [the original class](exact_problem.md). This page specifies a sampled-data procedure mathematically; sampled-data fitting remains disabled in this continuation. Only its deterministic weights and analytic fixtures were implemented and checked.

Split each source stratum by a fixed observation-index rule into a coefficient sample of size b_a and an independent evaluation sample of size m_a. Positive b_a for all five strata and positive m_0,m_3,m_4 suffice for the following construction. The e1/e2 evaluation observations can instead be assigned to coefficient learning in a predeclared allocation. No random split or sample generation is needed. Counts, delta and a failure probability determine all tuning below; target observations are never used.

**PROVED HERE — all four coefficients, using the verified P2/P7 construction.** Put `b_*=min b_a`, `epsilon_c=sqrt(log(10/a_c)/(2b_*))`. Use the historical separation constants

`a0=.5 exp(.5)`, `b0=2e`, `C=log(8)+.5`,

`R=sqrt((sqrt(1+4(C/delta)^2)-1)/2)`, `q=exp(-C)`,

`gamma=a0/2`, `J=b0(1+q)`.

Choose a finite telescope K_c and a coefficient grid with max covering distance tau in `[.5,1]^2`. Only finitely many profile queries are needed at the two calibration coordinates `z=-R,R` and the grid shifts. The private-anchor empirical quantiles yield these values by `sum_(k=1)^K_c Dhat_j(x-k)`. Take `T=R+K_c+2`, `U=T+1`. If `epsilon_c > Phi(-T)-Phi(-T-1)`, return zero for the target and a full-range risk bound. Otherwise set

`E_Q=epsilon_c * 2e[exp(U)+exp(h(U))h'(U)]/phi(U)`,

`P_c=4 K_c exp(R) E_Q + 2e(1+q) exp(-K_c)`, `O_c=exp(R)E_Q`,

`rho_c=min(.5, [2(P_c+O_c)+J tau]/gamma)`.

For each i=1,2, choose the lexicographically first grid minimizer of the maximum scaled residual at -R,R, with row scales `exp(R)` and `exp(-h(R))`. The five DKW bands hold simultaneously with probability at least `1-a_c`. The quantile bracket gives `E_Q`, telescoping gives P_c, and row dominance in P2 gives `||chat-c||_infinity<=rho_c`. Thus `E||chat-c||_infinity^2<=rho_c^2+a_c/4`. This is a finite, observation-only construction with explicit conditioning; it is expensive and conservative, and does not claim a fast coefficient rate. Whole profiles need not be reconstructed or stored for target evaluation. No optimization was executed.

**PROVED HERE — integrate weights before profiles.** Conditional on the coefficient split, fix K,B and use the combined A0,A3,A4 from [target_identity.md](target_identity.md), evaluated at chat. For `p=Phi(z)`, define

`w_a(p)=chi_B(z) A_a(z)/phi(z)`, `a=0,3,4`,

extended by zero outside `[Phi(-B-1),Phi(B+1)]`. The direct statistic is

`Jhat=sum_(a=0,3,4) sum_(i=1)^m_a Y_(a,(i)) integral_((i-1)/m_a)^(i/m_a) w_a(p) dp`.

These are order statistics of actual outcomes and weights computed from chat, delta and the calibrated design. The formula never requests true quantiles or profiles. Integrating a step empirical quantile gives this finite sum exactly. Its true-law counterpart is `J(chat,g)=sum_a integral w_a Q_a`. Because private-anchor differences do not involve the four coefficients, the identity remains valid when chat replaces c in the weights: `I(chat,g)=J(chat,g)+telescoping remainder+taper remainder`.

**STANDARD RESULT APPLIED — integrated empirical-quantile error.** Let `M_a=||w_a||_infinity` and `L_a=||w_a'||_infinity`. Write C_a for the source CDF and let `H_w(p)=integral_0^p w(u)du`. For positive outcomes,

`T_w(C)=integral_0^infinity [H_w(1)-H_w(C(y))]dy`.

Subtract its empirical version and apply Taylor's theorem:

`T_w(Chat)-T_w(C)=-integral w(C(y))(Chat-C)(y)dy+R`,

`|R| <= (L_a/2) integral (Chat-C)^2 dy`,

`E|R| <= L_a E|Y-Y'|/(4m_a) <= L_a mubar_a/(2m_a)`.

The expectation identity uses `E(Chat-C)^2=C(1-C)/m_a` and `integral C(1-C)=E|Y-Y'|/2`. This is a direct integrated quantile calculation, not the old uniform extreme-quantile bound. The weight is C1 with finite support in probability, and first moments suffice for these integrals.

For variance, replacing one observation changes the L-statistic by at most `M_a |Y_i-Y_i'|/m_a`: the difference is bounded by `M_a` times the empirical Wasserstein-1 distance, using the coupling that replaces only that observation. Efron–Stein therefore gives

`Var(T_w(Chat)) <= M_a^2 Var(Y)/m_a <= M_a^2 M2bar_a/m_a`.

Here fully known ceilings are

`mubar_a=2[exp(s1a)m_(1,1)+exp(s2a)m_(1,2)]`,

`M2bar_a=8[exp(2s1a)m_(2,1)+exp(2s2a)m_(2,2)]`,

with `(s10,s20)=(0,0)`, `(s13,s23)=(1,0)`, `(s14,s24)=(0,1)`. Gaussian moment ceilings from exact_problem.md can replace each m. Conditional independence across source strata now proves

`E[(Jhat-J)^2 | coefficient sample] <= sum_a M_a^2 M2bar_a/m_a + [sum_a L_a mubar_a/(2m_a)]^2`.

This is a finite-sample inequality for every profile in the original class and every fixed conditional chat,K,B. It is not an efficiency or sharp-rate assertion.

**PROVED HERE — computable ceilings for the weights.** Set `S=B+1`, `H1=1.1+.2S`, `A=8(1+H1)`, `D=4K(1+H1^2)+1.6`. The lattice bound in target_identity.md gives `|A_a|<=A`. Since `|f_j'|<=1`, `|h''|<=.2`, differentiation gives `|A_a'|<=D`. The first inequality for f2 follows from `|f2'(x)|<=phi(z)(|z|+.2)` with z=h_inverse(x). Consequently, uniformly over every chat in the coefficient box,

`M_a <= Mplus=A/phi(S)`,

`L_a <= Lplus=[1.5A+D+SA]/phi(S)^2`.

These explicit ceilings avoid an oracle supremum or an unproved numerical optimization of the constants. Tighter known-weight constants can be used only with adequate calculation error control.

**STANDARD RESULT APPLIED — reference covariance, explicitly.** The centered influence of one source L-statistic is

`ell_a(y)=-integral w_a(C_a(v))[1{y<=v}-C_a(v)]dv`.

Its covariance kernel in latent coordinates is `Phi(min(z,t))-Phi(z)Phi(t)`. Thus its variance is the double integral of that kernel against

`chi_B(z) A_a(z) F_a'(z)/phi(z)` and `chi_B(t) A_a(t) F_a'(t)/phi(t)`.

For reference stratum 0, use `A0=-A3-A4` inside this expression, retaining both cross terms. The true F' and exact influence variance are oracle quantities for analysis, not inputs to the estimator or its displayed envelope bound. The two anchor-reference differences are not independent. Independent coefficient splitting fixes their weights for this argument; no unjustified independence across reused observations is assumed.

**PROVED HERE — coefficient and nonlinear effects.** For every g in the original class, the ordinary target obeys the global bound

`|I(chat,g)-I(c,g)| <= C_c ||chat-c||_infinity`,

`C_c=4(e^2+e)(m_(1,1)+m_(1,2))`.

Proof: differentiate the two shifted terms involving each coefficient and use `|g'|<=2exp(x)` along the entire segment between c and chat. Sum the four derivative bounds and integrate the segment. This is a global finite-difference inequality, so no unknown second-order Taylor remainder is suppressed. Unknown profiles enter only through uniform envelopes; there is no additional fitted-profile radius in the evaluation bound.

Let `Rplus=10(e-1)^2(m_(1,1)+m_(1,2))exp(-K)`, `Tplus=C_T(delta)exp(-B)`,

`Vplus=sum_(a=0,3,4) Mplus^2 M2bar_a/m_a`, `bplus=sum_(a=0,3,4) Lplus mubar_a/(2m_a)`.

Clip Jhat to `[-Mbar,Mbar]`. In exact arithmetic, Minkowski's inequality and the preceding arguments give the complete unknown-coefficient bound

`MSE <= min{4Mbar^2, [sqrt(Vplus)+bplus+Rplus+Tplus+C_c sqrt(rho_c^2+a_c/4)]^2}`.

If a separately justified numerical target error is at most e_num, add e_num inside the square bracket. Every other displayed constant is known from counts, delta, the coefficient box and the profile envelopes. This proves feasibility and conservative finite risk, not the desired sharp order.

**PROVED HERE — source-only regularization and consistency.** A fully specified sufficient schedule is: let `n_*=min(b_a,m_0,m_3,m_4)`; choose `K_c=K=max(1,floor(.25 sqrt(log n_*)))`, `B=max(1,.25 sqrt(log n_*))`, `a_c=1/n_*`, and grid covering distance `tau<=1/n_*`. For small n where these definitions or the central-quantile condition fail, return zero with the full-range bound. No data-dependent target tuning is used. At fixed delta, the condition eventually holds. The coefficient stochastic term has the historical polynomial decay while its truncation term is O(exp(-K_c)). The direct-stage variance is `n_*^(-15/16+o(1))`, and its bias ceiling is `n_*^(-15/16+o(1))`, because Mplus squared and Lplus grow like polynomials times `exp((B+1)^2)`. The taper and telescope biases are O(exp(-.25 sqrt(log n_*))). Therefore

`MSE = O_delta(exp[-.5 sqrt(log n_*)])`.

Constants and onset are not uniform as delta tends to zero. This recovers the old sufficient slow order while avoiding target-stage whole-profile reconstruction. It does not improve the established original-class rate. An arbitrary data-selected K,B would require a uniform selection analysis; none is presumed here.

**PROVED HERE / OPEN — numerical feasibility versus present diagnostics.** Every sample coefficient weight is an integral of the known function `chi_B A_a` over a finite z interval. Its derivative has bound `1.5A+D`, so midpoint integration with mesh h has error at most `(1.5A+D)(2S)h/4` per interval. In a hypothetical execution, choose weight tolerances so `sum_(a,i) |Y_(a,(i))|*weight_error_(a,i)<=e_num`; this is possible with finite positive observations and uses no unknown labels. Known-warp bisection uses input-dependent brackets, avoiding the frozen pilot's fixed-bracket tail problem. Bracket, endpoint, exponential/CDF and arithmetic errors must also be enclosed, or retained separately. An approximate coefficient-grid minimizer with residual error at most xi_c adds `2xi_c/gamma` to the coefficient bound, along with any profile-query arithmetic errors.

The present implementation supplies deterministic float64 kernels, weights, taper and inverse, not a certified interval implementation of the sampled estimator. Its tests do not certify e_num. No sampled-data method was executed or enabled. The feasible statistical theorem is stated in exact arithmetic; a finite-precision guarantee needs the described error accounting rather than an assumption that solver tolerance is a probability guarantee.
