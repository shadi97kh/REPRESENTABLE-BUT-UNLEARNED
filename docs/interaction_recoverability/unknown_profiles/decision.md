# Unknown-profile extension decision

**incremental-known.** The completed construction establishes finite-action identification and a conservative nonparametric consistency bound under known warps, calibrated private anchors and explicit exponential profile envelopes. It is an ordinary inverse-problem plug-in estimator. No target-specific advantage or distinct learned ML contribution has survived this investigation, so no fitted pilot or training authorization is requested.

**Original rate: supported, allocation clarified — STANDARD RESULT APPLIED.** The old matching order remains
\(\min\{1,(n\delta^2)^{-1}\}\) for n observations **per stratum**, with constants uniform for known delta in (0,.1]. Both low-information and shrinking-alternative high-information lower bounds concern the complete experiment. For arbitrary allocations the informative count is \(\min(n_{e_1},n_{e_2})\); total sample count alone is insufficient. The old protected O(1/n) contrast is a known observed-mean cancellation, generally absent with unknown profiles. No historical proof or incremental-known decision was overwritten. See [rate_audit.md](rate_audit.md).

**Strongest completed extension result — PROVED HERE.** The declared class contains two genuinely unknown analytic profiles and the original four uncertain coefficients, while keeping warps and separation known. Differences between the private-anchor and reference quantile curves identify unit profile increments. A telescoping sum recovers each profile with an explicit lower-tail remainder. Two sufficiently separated latent quantiles then identify the coefficients through a derived diagonal-dominance bound. Therefore I12 is globally identified from the same five isolated action laws. This is not an inference from a zero score kernel.

The feasible reconstruction bound includes observed-quantile noise, estimated-profile error, profile truncation, coefficient discretization, numerical integration and outcome tails. For one explicit schedule it gives exact-arithmetic MSE
\(O_\delta(\exp[-\tfrac12\sqrt{\log n_*}])\) for fixed positive delta, \(n_*=\min_a n_a\). This is only a conservative consistency upper bound; it is neither minimax nor a uniform weak-separation rate. The old rate is a valid exponential-subclass lower bound, leaving a large unresolved gap. The full specification and proof are ordinary attachable Markdown: [model.md](model.md), [proofs.md](proofs.md), [estimator.md](estimator.md).

**Exact difference and remaining proof gap — OPEN.** DExtrI already learns profiles and warps, and its population theorem uses different profile/coverage restrictions. Our private-anchor difference identity handles an infinite profile class with only five action laws under stronger known calibration. This special-case distinction does not establish an ML advantage. The estimator reconstructs nuisance components and has no proven benefit over an equally informed inverse-problem control. Trabs' adjoint-range machinery already characterizes the relevant regularity question. The [novelty matrix](novelty_matrix.md) compares actual assumptions and targets, including bilinear completion and the two engression analyses.

The precise missing step is to solve or refute
\(DI_\eta[v]=\sum_a\pi_a E_a[\psi_aS_{\eta,a}(v)]\)
with finite weighted L2 norm for all admissible directions at fixed positive delta. The score and target derivative have been derived. There is no exact admissible null direction. A constructed sequence makes full-profile inversion unbounded in the specified weighted Sobolev norm, but its target derivative shrinks too; it does not establish failure of target regularity. A direct-target estimator would additionally need feasible nuisance estimation and nonlinear remainder control. No root-N claim or claim of impossibility at fixed delta is made.

**Actual new diagnostics — NUMERICALLY CHECKED.**

| Check | Output |
|---|---:|
| Unknown-profile telescoping identity, 36 exact-function evaluations | Maximum absolute residual 3.2862601528904634e-14 |
| Fixed-outcome log-density score derivative, 15 points | Maximum absolute difference 1.2624593037635634e-10 |
| Five score means integrated on [-12,12] | Maximum absolute value 6.46159083812943e-17 |
| Normalized score norms for profile-tail directions M=2,4,6 | .0610862, .0143758, .000656964 |
| Absolute target derivative / score norm for the same directions | .170527, .00797091, .000102337 |
| New correctness tests | 7 passed in .43 seconds |
| Samples drawn / sample fits / training runs / GPU use | 0 / 0 / 0 / 0 |

For K=2, B=2 and 95% simultaneous confidence, this proof's sufficient per-stratum sample threshold has log10 value 38.5545 at delta=.1 and 178.3293 at delta=.01. Both preflight bounds fail at n=1,000. These enormous thresholds describe this conservative procedure, not necessary statistical sample complexities. They rule out presenting this prototype as a useful demonstrated learner. Original four-delta instability outputs were reused, not rerun. All floating-point and quadrature errors are implementation diagnostics rather than certified mathematical enclosures.

The [diagnostic JSON](../../../runs/interaction_recoverability/prep-20260912-082154/unknown-profiles-v1/diagnostics/diagnostics.json) retains errors, intervals and conventions. The [command/test record](../../../runs/interaction_recoverability/prep-20260912-082154/unknown-profiles-v1/checks-v1.json) records actual commands. The [new implementation](../../../interaction_recoverability/unknown_profiles/core.py) has a default-disabled sample interface and no fit CLI. No sample estimator was enabled, and its end-to-end fitted behavior is untested. [Correctness tests](../../../tests/interaction_recoverability_unknown_profiles/test_math.py) use only exact-function fixtures and analytic bounds.

**Preservation and resources.** The [verification record](../../../runs/interaction_recoverability/prep-20260912-082154/unknown-profiles-v1/verification.json) checks both the original preservation manifest and the continuation's snapshot of prior recoverability artifacts. Existing code, experiments, negative results, user files and the old pilot are retained. No AGENTS.md was found in the applicable searched locations. The [resource notes](resources.md) explain charges, unmeasured discovery overhead and the retained clean-stop reserve. The authoritative [ledger](../../../runs/interaction_recoverability/prep-20260912-082154/ledger.jsonl) was extended without reset; [snapshot 26](../../../runs/interaction_recoverability/prep-20260912-082154/resource_snapshot-26.json) includes the final verification job. CPU, aggregate local job wall time and conservative charges are separate; elapsed research/conversation time is not local compute-job wall time. One CPU, one numerical thread and zero GPU were retained.

The [estimator document](estimator.md) gives reproducible preparation commands. No fitted comparison is proposed because no credible constructive advantage survived. The original known-profile pilot remains frozen and unexecuted. The [model document](model.md) states the measured siRNA contrasts, five-group source/anchor requirements, held-out joint evaluation and distinction from frozen-GNN explanations; synthetic fixtures supply no biological validation.

**One next action:** resolve the fixed-delta I12 adjoint-range problem in P6 against an equally informed direct-target inverse-problem baseline, before requesting any fitting budget.
