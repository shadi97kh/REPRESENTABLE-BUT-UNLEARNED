# Learning a selected interaction from a finite action dictionary

The primary target, fixed before the diagnostic was executed, is

\[
I_{12}=\mu(e_1+e_2)-\mu(e_1)-\mu(e_2)+\mu(0).
\]

The secondary candidate contrast is \(\Delta=\mu(a_\star)-\mu(b_\star)\); when needed below, \(a_\star=e_1+e_2\) and \(b_\star=0\). These are means on the original positive response scale, in arbitrary synthetic units. A logarithmic transformation would define a different interaction and is not performed. Background and assay context are fixed; actions are four declared coordinates in \([0,1]^4\), with a finite observed dictionary \(\mathcal A_{\rm obs}=\{0,e_1,e_2,e_3,e_4\}\).

**PROVED HERE — observation model.** For each observed action \(a\), the complete experiment contains \(n_a\) independent outcomes \(Y_{a,i}=G(a,Z_{a,i})\), where all \(Z_{a,i}\sim N(0,1)\) are independent across observations and actions. The latent values are unobserved. A single outcome uses the same latent value in its two summands. No cross-action latent pairing is observed. Thus the experiment law is \(\mathbb P_G=\bigotimes_{a\in\mathcal A_{\rm obs}}P_{G,a}^{\otimes n_a}\). Action counts are fixed in advance; adaptive acquisition would require a different proof. Numerical displays use \(n_a=32\) solely when evaluating analytic bounds; no sample is generated.

**PROVED HERE — normalized class.** Fix a *known* \(\delta\in(0,0.1]\), known Gaussian latent law, known profiles \(g_1=g_2=\exp\), and known warps \(h_1(z)=z\), \(h_2(z)=z+\delta z\sqrt{1+z^2}\). Define

\[
G_c(a,z)=e^{a_3+z}\,u_1^{a_1}u_2^{a_2}
 +e^{a_4+h_2(z)}\,v_1^{a_1}v_2^{a_2},
\qquad (u_1,u_2,v_1,v_2)\in[\sqrt e,e]^4.
\]

Only these four coefficients are unknown. The supports are \(\{1,2,3\}\) and \(\{1,2,4\}\). Fixed unit loadings on the private anchors \(a_3,a_4\) label components and remove permutation and index-rescaling ambiguities. The noise median is zero, profiles and additive term are fixed, and there is no unfixed additive-constant exchange. Parameter error concerns these normalized coefficients, not an arbitrary neural parameterization. The supplied T and Q models are two members of this class.

Write \(m_0=e^{1/2}\) and \(m_\delta=E e^{h_2(Z)}\). They are known model-class integrals, not fitted biological quantities. The primary target is

\[
I_{12}(c)=m_0(u_1-1)(u_2-1)+m_\delta(v_1-1)(v_2-1).
\]

The statistical question is the finite-sample minimax MSE for this functional. Computational implementation of a two-by-two inverse is a separate question. Full recovery of all coefficients is not the definition of success.

**STANDARD RESULT APPLIED — completed restricted rate.** For balanced counts \(n_a=n\ge1\), the proof notes establish constants \(0<c<C<\infty\), independent of \(n\) and \(\delta\in(0,.1]\), such that

\[
c\min\{1,(n\delta^2)^{-1}\}\le
\inf_{\widehat I}\sup_{c\in[\sqrt e,e]^4}E_c(\widehat I-I_{12}(c))^2
\le C\min\{1,(n\delta^2)^{-1}\}.
\]

This uses a complete-experiment likelihood lower bound and a feasible empirical-quantile inverse. It is a concrete weak-identification calculation with known profiles, rather than a new learning mechanism. Constants are conservative, and no optimal-constant claim is made. The parameter symbol \(c\) in the supremum denotes the coefficient vector; the two unnamed constants in the display are unrelated to it.

**PROVED HERE — protected target.** For the different, secondary contrast

\[
J=\mu(e_1+e_3+e_4)-\mu(e_1)-\mu(e_3+e_4)+\mu(0),
\]

the class identity \(J=(e-1)[\mu(e_1)-\mu(0)]\) permits ordinary sample-mean recovery at MSE \(O(1/n)\), uniformly as \(\delta\downarrow0\). Neither component needs to be recovered. This is an explicit example of an estimable functional avoiding weak directions; it remains secondary and does not replace the hard primary target after seeing results.

**OPEN — intended broader contribution.** Unknown warp functions, unknown profiles, unknown separation, different latent dimensions, finite categorical chemistry without calibrated anchors, and misspecified interaction structure are outside this result. A theorem using a uniform inverse inequality as an assumption would merely encode the desired answer. The missing distinction is an observable, non-tautological condition and an attainable target bound when those nuisance mechanisms must themselves be learned.

The source population theorem has a different observation restriction. DExtrI uses conditional laws across intervals on each axis; a few categorical actions do not supply that coverage. In this restricted class, two quantiles at each of two isolated observed actions suffice only because the latent basis is already known. No continuous-axis theorem is silently transferred to binary chemistry. See the [novelty comparison](novelty_matrix.md) and [proofs](proofs.md).
