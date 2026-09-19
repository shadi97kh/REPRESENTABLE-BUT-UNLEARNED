# Prior overlap and the remaining contribution question

This is a proof-only, bounded comparison against the requested primary results. It is not an exhaustive priority search or external independent review. The mathematical status in [decision.md](decision.md) governs the construction; this comparison cannot supply a missing attainment proof.

## Exact object being compared

**PROVED HERE in the accompanying proofs:** the additional model-specific result concerns five isolated source laws, a known Gaussian latent distribution, known private unit anchors and a known warp
\[
h_{\delta,\theta}(x)=x+\delta x\sqrt{1+\theta x^2},
\qquad \theta\in[1/2,2],\quad 0<\delta\le 1/3200000.
\]
The profiles retain the original analytic derivative envelope, rather than a neighborhood whose radius shrinks with \(\delta\). All four coefficients remain unknown in \([.5,1]\). The conclusion is source-data pointwise attainment of the interaction at each fixed truth and fixed known warp, with an explicit, generally poor constructed-influence bound. See [problem.md](problem.md), [commutator_and_calibration.md](commutator_and_calibration.md), [attainment.md](attainment.md) and [conditioning.md](conditioning.md) for the actual quantifiers and proof obligations.

The interval is a sufficient small-parameter interval. It does not contain the historical \(\delta=.1\) example and is not a necessary restriction on statistical recovery. Each known pair \((\delta,\theta)\) indexes a separate experiment. No conclusion here estimates an unknown warp or permits a sequence \(\delta_N\to0\) without additional analysis.

## Requested primary results

### Trabs: adjoint ranges and information bounds

**STANDARD RESULT APPLIED.** [Trabs, arXiv:1307.6610v2, Theorem 2.7](https://arxiv.org/html/1307.6610v2) characterizes regularity in locally regular indirect experiments through the adjoint-score range. The adjoint pseudoinverse gives the efficient influence and convolution bound. Example 2.9 applies the framework to nonlinear forward operators; the following paragraph does not supply general nonlinear attainment.

Here, deriving an influence against the complete five-law score operator verifies that established framework. Additional work consists of observable inversion, calibration with four unknown coefficients, and a source-data statistic whose growing-series and nonlinear remainders are negligible. A nonminimum-norm influence does not establish efficiency. The scope of Trabs's nonlinear discussion is not evidence of priority.

### Kaji: inverse quantiles and integrated statistics

**STANDARD RESULT APPLIED.** [Kaji, arXiv:1910.07572v1, Theorem 4.1 and Proposition 4.2](https://arxiv.org/html/1910.07572v1) establish integrable inverse-quantile differentiability and Gaussian limits under positive-density and integrated-tail conditions. Theorem 5.1 treats an integrating-function class with a fixed Lipschitz bound.

The present proof must verify the warp-dependent moments and control the growing contraction depth, coefficient-dependent queries, transport joins and changing tail cutoff. A fixed integrated-quantile limit alone does not control those operations. The resulting remainder argument is additional verification using established empirical-process machinery, not a new quantile limit principle.

### DExtrI: structural identification from interval information

**STANDARD RESULT APPLIED AS A COMPARISON.** [DExtrI, arXiv:2608.19849v1, Theorems 3.3–3.4](https://arxiv.org/html/2608.19849v1) identify and extrapolate warped ridge models using conditional laws across coordinate-axis intervals under A1–A7. Theorem D.3 covers broader profiles conditional on rigidity and other hypotheses. Theorem D.9 proves rigidity under a complex-singularity assumption.

Exactly five isolated laws with known warp and private anchors supply different information. Population identification under interval observations does not itself give the present sampling theorem. Conversely, leaving A7's profile class is not proof of novelty: D.3 requires separate consideration.

Condition S in D.9 concerns nonempty discrete nonremovable complex singularities of profile derivatives, with conditions on rescaled profile differences. Branch points of the square-root warp do not establish S. The original profile envelope does not imply S. Nor should distinct known \((\delta,\theta)\) values be pooled into an unknown-warp comparison to manufacture a failure of DExtrI's growth-order assumptions.

### Bennett and collaborators: functionals under weak nuisance identification

**STANDARD RESULT APPLIED AS A COMPARISON.** [Bennett et al., arXiv:2208.08291v3, Theorems 2–3 and Section 4.1](https://arxiv.org/html/2208.08291v3) treat functional recovery despite weak or absent identification of the full conditional-moment nuisance. Their strong identification condition is \(\alpha\in\operatorname{Range}(P^*P)\); existence of a square-integrable debiasing nuisance uses the weaker adjoint-range condition. Theorem 3 provides a doubly robust bias identity.

The nonlinear five-law response experiment is not automatically their conditional-moment model. Observable inversion and attainable remainder control must still be demonstrated. However, recovery of a target without accurate recovery of every nuisance is already an established principle and cannot be presented as this project's new contribution.

### Yao and de la Llave: cohomological equations and series

**STANDARD RESULT APPLIED.** [Yao and de la Llave, arXiv:2110.15893v1, equations (15)–(18) and Section 3.4](https://arxiv.org/html/2110.15893v1) discuss cohomological equations \(\phi=l\,\phi\circ a+\eta\), weighted composition-series solutions, convergence conditions and accelerated summation. Algorithm 2 repeatedly composes and combines partial sums. Their discussion also recognizes rounding and truncation limitations.

The centered profile equation has multiplier one, so its normalization and function space matter. After subtracting the common fixed-point value, composition by \(C\) contracts a Hölder seminorm by \(\lambda^\gamma\), giving the usual Neumann bound \((1-\lambda^\gamma)^{-1}\). Alternatively differentiation introduces multiplier \(C'\). Thus an unrestricted multiplier-contraction theorem cannot simply be cited with \(l=1\), but the normalized inverse is established contraction machinery. The quantile kernel's one-half-Hölder increment bound gives \((1-\sqrt\lambda)^{-1}\). Neither the series nor this geometric denominator is a new principle.

## Exact additional statistical work

**PROVED HERE:** two changes address the old construction's dependence on a specially calibrated neighborhood.

1. The six-term observable commutator supplies central inversion across an explicit family, with fixed compact central queries and a derived contraction margin proportional to \(\delta\).
2. Calibration at known \(\pm R\), where \(h(R)-R=3\), gives a uniform population separation argument over the original profile envelope. This avoids shrinking the profile class with \(\delta\), while introducing an explicitly charged extreme-quantile cost.

The construction still uses only the existing source laws. Distinct rows of an observable population calibration map, diagonal dominance, projected minimum-distance iteration, quantile inversion, telescoping, and target integration remain established tools. The additional model-specific result is their statistical combination, proved in the accompanying documents, in this experiment, including simultaneous schedules and all five influence components.

Every equally informed inverse-functional, direct-target, likelihood or confidence-set method may use the new identity and calibration. A control must not be forced to reconstruct complete profiles accurately when its target does not require that. This continuation proves no comparative superiority over such controls.

## Conditioning is not a matched separation rate

**PROVED HERE, only with the scope proved in [conditioning.md](conditioning.md):** a sufficient bound for the constructed influence can grow exponentially in inverse separation. The contraction-series factor is polynomial, but the calibration location obeys \(R^2\asymp\delta^{-1}\), so Gaussian quantile norms introduce an exponential factor. Pointwise root-\(N\) attainability remains consistent with very unfavorable constants.

The accompanying complete-experiment exponential-subclass lower argument is rederived for this warp family and its allowed parameters. A verified squared-risk lower order \(\min\{1,(N\delta^2)^{-1}\}\) is a worst-case statement for that subfamily; it does not match an exponential sufficient upper constant. The old numerical constants and the historical warp restriction cannot simply be reused. A large upper bound for this construction is not an impossibility theorem, and a two-point lower bound is not an attaining estimator.

**OPEN:** sharp separation dependence, useful triangular-array guarantees as \(\delta_N\to0\), efficient variance, finite-sample accuracy, mean-squared-error attainment and superiority over equally informed controls. No unmatched order is described as optimal.

## Significance and preservation

**OPEN:** publication priority and substantive significance. Extending one calibrated point to a stated warp family and the original profile envelope is the concrete model-specific strengthening proved in the accompanying documents. It does not alone establish a new general ML principle. Sharpness, a meaningful uniform separation regime, or a demonstrated necessity or benefit of the observation design would strengthen the contribution case. Their absence is a significance limitation, not an automatic contradiction of a pointwise theorem.

The failed learned-target pilot, unfavorable historical influence scale, backend no-go, unresolved efficiency decision and absence of biological validation remain unchanged. No neural training is needed merely to assess a theorem, and no biological siRNA claim follows from these five-law assumptions.

## Source-inspection record

Cached primary text inspected:

- `runs/interaction_recoverability/prep-20260912-082154/unknown-profiles-v1/trabs_v2.txt`, lines 315–390 and 426–474.
- `runs/interaction_recoverability/prep-20260912-082154/sources/dextri_v1.txt`, lines 329–512, 1083–1128, 1520–1585 and 1924–1984.
- `runs/interaction_recoverability/prep-20260912-082154/sources/strong_functionals.txt`, lines 575–745.

Targeted primary HTML inspected: Kaji v1 Sections 4–5; Yao and de la Llave v1 Section 2.2.1.2, equations (15)–(18), Section 3.4, equation (22), Algorithm 2 and Remarks 16–21. No requested source was unavailable. We did not independently reprove every cited external theorem or conduct an exhaustive search. The comparison used read-only inspection and source retrieval; no numerical or project code was executed and no ledger was invoked.
