> Publication copy of `docs/interaction_recoverability/unknown_profiles/model.md`. Mathematical and result content is retained; execution material is omitted. This dated record retains its original conclusions and limitations, including conclusions superseded by later stages. Original SHA-256: `bb4e9e43e7b357df20546118760ddceb4eed877c8d39d12eac091c53967df36c`.

# Unknown profiles, fixed known warps: model v1

This is a new class, not a change to the historical known-profile result. All structural restrictions below are assumptions; they are not learned from siRNA data.

**PROVED HERE — declared experiment and normalization.** Fix known \(0<\delta\le .1\), \(h(z)=z+\delta z\sqrt{1+z^2}\), and independent unobserved \(Z_{ai}\sim N(0,1)\). For \(a\in[0,1]^4\), set

\[
F_a(z)=g_1(z+\alpha_1a_1+\alpha_2a_2+a_3)
       +g_2(h(z)+\beta_1a_1+\beta_2a_2+a_4),
\quad (\alpha_1,\alpha_2,\beta_1,\beta_2)\in[1/2,1]^4.
\]

Both profiles are unknown members of the following infinite-dimensional class:

\[
\mathcal G=\{g:\mathbb R\to(0,\infty):g\text{ real analytic},\quad
g(0)=1,\quad g(-\infty)=0,\quad
\tfrac12e^x\le g'(x)\le2e^x,\quad
|g''(x)|,|g'''(x)|\le10e^x\}.
\]

There are no additive terms, unknown intercepts, unknown warps, or unknown separation parameters. Supports remain \(\{1,2,3\}\), \(\{1,2,4\}\). The private loadings are known to equal one and label the components. Known warps fix index scale and location; the boundary and value normalization remove constant exchange. These are substantial restrictions, not cost-free biological calibration.

The complete observed experiment is still
\(\bigotimes_{a\in\{0,e_1,e_2,e_3,e_4\}}P_{\eta,a}^{\otimes n_a}\), with fixed counts and independent outcomes across all strata. A response's two summands share its one latent draw; there is no observed pairing across actions. Let \(N=\sum_a n_a\). Balanced \(n\) means **per stratum**, so \(N=5n\). The constructive bound uses \(n_* =\min_a n_a\); positive counts at all five strata are its sufficient design. Population identification refers to knowledge of these five whole laws, not finitely many observed outcomes or laws over continuous action intervals.

**PROVED HERE — genuinely unknown functions and moments.** Integrating the derivative bounds gives \(e^x/2\le g(x)\le2e^x\). The class includes \(e^x(1+\epsilon\sin(\omega x))\) for every \(\omega\in[1,2]\), \(|\epsilon|\le .05\), as well as sufficiently small analytic perturbations with bounded weighted derivatives. The distinct frequencies already give infinitely many linearly independent possibilities; the estimator below uses no frequency dictionary. Neither function is fixed to its generating value. The exponential subclass is retained exactly via \(u_i=e^{\alpha_i},v_i=e^{\beta_i}\).

The derivative lower bounds make every \(F_a\) strictly increasing onto \((0,\infty)\). Thus \(F_a(z)=Q_a(\Phi(z))\), and the density is \(p_a(y)=\phi(F_a^{-1}(y))/F_a'(F_a^{-1}(y))\). Comparability to exponentials gives all response moments of order \(p\) with \(p\delta<1/2\), uniformly for \(p=1,2\) and allowed delta. It does not give bounded outcomes or a fifth moment at \(\delta=.1\).

**PROVED HERE — fixed primary target.** On the original positive response scale,

\[
I_{12}=E\{F_{e_1+e_2}(Z)-F_{e_1}(Z)-F_{e_2}(Z)+F_0(Z)\}.
\]

For each component the integrand is
\(g(x+c_1+c_2)-g(x+c_1)-g(x+c_2)+g(x)\), with \(x=z\) or \(h(z)\). It need not be positive because convexity was not assumed. If \(m_0=e^{1/2}\), \(m_\delta=E e^{h(Z)}\), then
\(|I_{12}|\le 2(e+1)^2(m_0+m_\delta)\). These known exponential moments bound unknown profiles; using them is not substituting oracle profiles. Loss is squared error for this target. Recovering all profiles on the whole line is a stronger task and is not the success criterion.

**PROVED HERE — global identification.** The explicit reconstruction and coefficient separation in [proofs.md](proofs.md) identify the entire normalized model, and hence this contrast, from these five laws for every fixed positive delta. No indistinguishable pair exists within this class. A score calculation alone is not used as proof of this fact. At delta zero the retained exponential subclass still has the old nonidentification example; no theorem includes that boundary.

**OPEN — application assumptions.** For siRNA, specify two particular modifications at two positions on the same sequence/backbone, endpoint, dose, time, cell system and assay batch. Biological \(I_{12}\) requires means for untreated/reference, each single change and their joint change, with independent experimental replicates and batch/randomization information. A held-out joint group can evaluate extrapolation but must be unavailable to every learner. This model additionally needs two measured private-anchor interventions whose unit loadings and known latent warps are scientifically justified; ordinary modification labels do not provide those assumptions. Without such calibration the theorem is not a biological guarantee. Replicate counts at all five source groups must be reported, and joint-group evaluation uncertainty must be retained. Frozen-GNN query differences assess that model's behavior, not measured biological interactions. No graph encoder, therapeutic claim or synthetic biological validation is supplied here.
