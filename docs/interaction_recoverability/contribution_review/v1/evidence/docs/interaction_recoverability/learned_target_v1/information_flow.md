# Information flow and access audit

Status: implemented, exercised by the development fit and frozen for the pilot.

| Object | Learner access | Use and restriction |
|---|---|---|
| Five action-labeled outcome arrays, counts | Every applicable method | The source loader rejects extra archive keys, latent arrays and metadata |
| Known Gaussian law, h_delta, delta, anchors, coefficient box, profile restrictions | Every applicable method | Shared structural information, not estimated or biologically validated |
| Target contrast formula | Every applicable method | Its numerical true value is withheld |
| Public integration nodes | Every applicable method | Integrate candidate models; never supplied as observed individual latent coordinates |
| Public optimization/split seed | Every applicable method | Independent of the private data-generation seed |
| Fit part F | Fit nuisance and geometry | No use of that rotation's E outcomes |
| Selection part C | Checkpoints, lambda, enrichment/finite paths | All choices frozen before E correction |
| Evaluation part E | Compute that rotation's final correction only | Cannot select its nuisance, direction space or lambda |
| All source parts | Pooled plug-ins and direct-target set control | Stronger pooled controls get the same total source observations |
| True profiles/shifts, family label, generator seed, true target | Evaluation module only | Separate restricted files; not passed to learner functions |
| Joint-action outcomes | No learner | Not generated as training data |

`experiment.py` imports no evaluation module. Its source interface accepts only an archive containing `y0,...,y4`, delta, a public initialization seed, frozen settings, and the common metric. AST/interface checks supplement inspection. The command dispatcher imports the evaluator only for development generation/scoring and the designated final generation/evaluation/plot stages. The generator never stores individual latent arrays. Private files have owner-only permissions; this is code/interface separation within one Unix account, **not an adversarial OS security boundary**.

The development dataset is independent of the final datasets. The only development statistic used to choose full versus minimum pilot is measured end-to-end runtime. Final config and source-code hashes are frozen before final sample generation. Private generation specifications have a public commitment hash. Final evaluation first hashes all saved predictions, then reads the target files. The training CLI refuses to run after the final evaluation directory exists. No post-evaluation parameter changes or refits are authorized by a favorable or unfavorable result.

The three cyclic estimates are dependent: an observation can be fit data in another rotation. Nevertheless each rotation's E outcomes are independent of its own F,C and corresponding correction. Conditional risk identities apply rotation by rotation. Jensen's inequality bounds the mean estimator's unconditional MSE by the mean of rotation MSEs without claiming independence. Per-fold variance diagnostics are not averaged into a nominal confidence interval.

Reconstruction of the full latent model is not the pilot metric. True-profile errors, known latent values, or oracle source-law distances do not select any fitted method. The target is the unobserved action contrast; the separate evaluator uses truth only after predictions are saved. Model-specific truth generation and the structural assumptions remain synthetic and do not validate siRNA biology.
