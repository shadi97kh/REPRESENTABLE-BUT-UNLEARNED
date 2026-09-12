> Publication copy of `docs/interaction_recoverability/estimator.md`. Mathematical and result content is retained; execution material is omitted. This dated record retains its original conclusions and limitations, including conclusions superseded by later stages. Original SHA-256: `9a7a38de07f9e3cbaa6435c0fe42b90d6b895ac4a071f67905a119d51d1ebad4`.

# Implementable recovery and reduction to known methods

**STANDARD RESULT APPLIED.** The implemented estimator is empirical-quantile minimum distance with a known two-column latent basis. It is a control for the normalized class, not a new architecture. The matched-rate proof is in [proofs.md](proofs.md); the complete observation model is in [problem.md](problem.md).

| Quantity | Information status |
|---|---|
| Background, original response scale, action labels and sample counts | Known design |
| Gaussian latent marginal, exponential profiles, private-anchor coefficients | Known structural assumptions, supplied equally to all applicable controls |
| Delta and its warp function | Known in this restricted theorem and proposed control pilot; oracle-only relative to an uncalibrated biological problem |
| \(m_0,m_\delta\), central quantile design and separation | Calculable from those known assumptions; numerical integration errors must be retained |
| Four amplitude coefficients and two observed-action quantile vectors | Unknown and estimated from outcomes |
| Latent Z for an individual outcome | Unobserved; never supplied to a feasible learner |
| Joint-action labels, generating coefficients and exact target | Evaluation-only oracle information |
| Full unknown profiles/warps or separation learned from these observations | OPEN; not covered by the estimator |

The implementable algorithm is:

```text
Input: independent outcomes at e1 and e2; known delta/profiles; fixed .25/.75 levels.
1. Take empirical inverse-CDF quantiles at the two fixed levels in each stratum.
2. Construct the known 2-by-2 matrix V_delta.
3. Solve V_delta c_i = empirical_quantiles_i for each i=1,2.
4. Clip each coefficient to [sqrt(e), e], the declared parameter box.
5. Return m0*(u1-1)*(u2-1) + m_delta*(v1-1)*(v2-1).
```

`estimate_from_samples` is implemented in [estimator.py](../../../../interaction_recoverability/estimator.py) and rejects calls unless `fitting_authorized=True`. No such call was made. Correctness checks supply exact analytic quantiles to the separate algebra interface; they neither sample outcomes nor fit a toy estimator. The protected contrast J uses only two sample means in a future application; its deterministic identity is checked without sample fitting.

Sorting costs \(O(n_{e_1}\log n_{e_1}+n_{e_2}\log n_{e_2})\) time and linear memory with the simple implementation; order-statistic selection could reduce time. Matrix construction/inversion and target evaluation are constant-dimensional. Moment quadrature is shared by all methods that use the same structural model and must be charged once plus reuse cost. No latent-model training or GNN is necessary.

**PROVED HERE — high-probability error statement.** With two balanced strata, \(\epsilon=\sqrt{\log(4/\alpha)/(2n)}\), the two DKW events hold simultaneously with probability at least \(1-\alpha\). When \(\epsilon\le1/8\), all coefficient errors are at most \(\min\{e-\sqrt e,L_\delta\|V_\delta^{-1}\|_\infty\epsilon\}\). Multiplying by \(K_\delta\) from the proof notes bounds target error. For smaller n the code declares this central-range bound uninformative; it does not relax the condition or substitute a fitted condition number. A confidence interval formed from this deterministic bound would be valid under the stated class but could be very wide. No useful-width or adaptive-coverage claim is made.

The main reductions are consequential:

| Candidate ingredient | Reduction and implication |
|---|---|
| Two-quantile recovery | Ordinary minimum-distance inversion; smallest singular value is of order delta |
| Likelihood competitor | \(p_c(y\mid a)=\phi(G_c^{-1}(y))/\partial_zG_c(G_c^{-1}(y))\); optimizing its product likelihood is ordinary MLE, with inverse/Jacobian costs included |
| Weighted quantile equations | Generalized least squares / efficient minimum distance once their covariance is known or consistently estimated; unknown covariance is not free oracle information |
| Selected contrast projection | For a known linear observation operator A, a linear target \(\ell^\top c\) is identifiable exactly when \(\ell\) lies in its row space; the implemented pseudoinverse check is this standard criterion |
| Protected target J | An exact observed-mean functional, so weak component directions cancel before inversion |
| Primary target I12 | A bilinear function of the two coefficient rows. A structural bilinear quantile-completion control performs the identical inverse and must match the proposed estimator |
| Rank-two completion using means alone | The baseline means on a cross need not identify two latent columns. Distributional observations are additional information; this does not prove that a particular new learner is required |
| Confidence set or profile interval | Known inference machinery; shrinking width without coverage is not an improvement |

The direct functional collision is [Bennett et al., arXiv 2208.08291v3](https://arxiv.org/abs/2208.08291v3): strong functionals of ill-posed or nonunique nuisance functions, operator range conditions and debiased estimation already have a developed theory. The present finite-dimensional calculation does not establish a weaker or sharper result. [Donoho's optimal-recovery analysis](https://web.stanford.edu/dept/statistics/cgi-bin/donoho/wp-content/uploads/2018/08/SEOR.pdf) likewise makes the difficulty of a target functional, rather than whole-object error, central.

**OPEN.** Estimating an unknown delta would alter both the inverse matrix and the target moment. Replacing it with its true value or a training-independent oracle estimate cannot support an adaptive claim. Unknown profiles would turn this into a nonlinear inverse problem requiring a new nuisance class, calibration design and error propagation argument. No guarantee here covers that case. In particular, the desired target modulus cannot simply be assumed as the quantitative condition.

The proposed unexecuted comparison is frozen in [the pilot configuration](../../../../configs/interaction_recoverability_pilot.json). Its purpose is to test a precise future claim if one is developed. Under the present known-warp result, equivalence to its bilinear/minimum-distance control is already decisive, so training authorization is not requested.
