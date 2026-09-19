# Paper tables

Generated 2026-09-05 by `scmp.export_results` from immutable result files. No number here was computed at export time; each is read from a file under `runs/`, and any file with a `.sha256` sidecar is re-hashed and refused if it moved.

Three kinds of result are kept apart. **A** are supporting lemmas, proved and checked by exhaustive enumeration. **B** are engineering checks. **C** are trained empirical results. A result in A or B is not evidence for a claim in C.


## A. Supporting lemmas (no trained parameters)

| lemma | measurement | value | source |
|---|---|---|---|
| Endpoint exactness iff aggregator monotone under multiset inclusion (either direction) | 7 aggregators x 25 seeded trials | holds: True | c3_aggregator_characterization.json |
| sum+ReLU under non-negative parameters is exactly affine | additivity violation | 2.157e-16 | c4_affinity.json |
| Relaxation gap equals a walk-count ratio kappa_L | max relative error vs exhaustive enumeration | 2.121e-16 | c7_kappa_characterization.json |
| Relaxation-gap exponent is set by network depth, not by the data | exponent per layer, depths 1-4 | 0.912-0.912 | c5_slack_scaling.json |
| Gap is unbounded on downward-closed, non-union-closed families | slack at m=3 -> m=48 | 9.0 -> 2290.3 | c5_slack_scaling.json |
| Gap is exactly 1.0 on join-closed families | brute-force verified | True | c5_slack_scaling.json |


## B. Engineering checks (no trained parameters)

| check | scope | result | source |
|---|---|---|---|
| Bound soundness, every intermediate bracket | 1,262,952 checks over exhaustively enumerated families | 0 violations | audit_bounds |
| Conditional-support soundness | 500 randomised cases | 0 violations | audit_conditional |
| Oracle correctness | counts vs Motzkin numbers, n=0..10 | exact match | audit_oracles |
| Refinement certificates | 72 certificates, independent checker | 72/72 accepted, 2222 obligations | check_certificates |
| Model-effect interval containment | 75 intervals vs enumeration | 0 failures | effects-20260905-161546.json |
| Reproduction of the winning development configuration | 35/35 source hashes, 5 data hashes | reproduced, 5/5 metrics | reproduction-20260905-172345.json |


## C. Trained empirical results (fitted parameters involved)

These are the only rows that bear on a predictive or head-to-head claim.

### Preregistered gates

| gate | point | 95% CI | threshold | status |
|---|---|---|---|---|
| G1 conditional vs unconditional | - | - | - | UNDETERMINED |
| G2 learned vs deterministic allocation | - | - | - | UNDETERMINED |
| G3 predictor noninferiority | -0.3361 | [-0.4858, -0.1869] | 0.05 | PASS |
| G4 residual predictive value | -0.2680 | [-0.4076, -0.1284] | - | FAIL |

Hardware: NVIDIA TITAN RTX (device 1, 0.15 memory cap). Budget used 0.2389 core-hours, within cap True.


### Verifier arms at matched wall time

| arm | n | median gap | mean gap | closed | timeout | median wall s |
|---|---|---|---|---|---|---|
| box | 330 | 0.0000 | 0.0008 | 99.1% | 0.9% | 0.060 |
| fixed_conditional | 330 | 0.0000 | 0.0135 | 83.6% | 22.1% | 0.567 |
| intermediate_unconditional | 330 | 0.0000 | 0.0008 | 99.1% | 0.9% | 0.060 |
| learned_conditional | 330 | 0.0000 | 0.0160 | 74.8% | 33.9% | 0.632 |
| terminal_only | 330 | 0.0000 | 0.1743 | 69.7% | 30.3% | 0.089 |

### Predictor arms, gene-disjoint development split

| arm | seeds | mean Spearman | sd |
|---|---|---|---|
| certifiable_gamma0 | 5 | +0.7253 | 0.0109 |
| certifiable_gamma1 | 5 | +0.4573 | 0.1593 |
| unconstrained | 5 | +0.1212 | 0.0653 |

Label semantics: the target is marginal over latent structures, so these rows measure a marginal-label fit only (`preregistration.md` section 6).



## Limits of validity

| dimension | what was actually covered | what is not claimed |
|---|---|---|
| Families | non-crossing pairings, layered-DAG paths, path matchings; ground sets up to 45 edges | no family whose membership needs more than the implemented oracles; no unbounded families |
| Architecture | 2-layer message passing, sum aggregation, ReLU, signed weights, linear readout, non-negative additive anchor | no attention, no normalisation, no gating, no depth beyond 5 (tested), no aggregation outside the monotone class |
| Numerics | float64 with no directed rounding; validated interval arithmetic exists and is tested but is NOT wired into the bound propagation | no bound in this work is numerically certified; all are `float_unverified` diagnostics (`docs/bounds_proof.md`) |
| Scale | exhaustively enumerable families for every soundness claim; windows up to 150 nt for descriptive statistics only | no soundness claim is validated beyond enumerable size |
| Statistics | 5 seeds, cluster bootstrap by instance | predictor CIs rest on 5 seeds and one held-out gene; they are wide and are not a population claim |


## Contribution, stated without the application domain

The object is a certified upper bound on a neural network's output over a
**combinatorially constrained family of input graphs**, where membership in the
family is decided by an exact oracle rather than by a norm ball or an edge
budget. Three parts:

1. **An additive anchor plus a bounded residual.** The model is
   `f(X, z) = b(X) + a(X)^T z + gamma * g(X, B + Pz)`. The anchor is exactly
   optimisable over the family by the oracle; only the residual needs relaxing.
   At `gamma = 0` the bound is exact by construction.
2. **Conditional message support.** A bilinear message term `z_e * h` is bounded
   using `h`'s maximum over the family *conditioned on* `z_e = 1`, which is
   never worse than the unconditional envelope and is strictly better whenever
   conditioning excludes competing structure.
3. **A characterisation of when relaxation is loose.** The gap equals a
   walk-count ratio; its exponent is set by network depth; it is unbounded on
   downward-closed families that are not union-closed and exactly 1 on
   join-closed ones.

Part 3 is the general result and is independent of any application. Parts 1 and
2 are mechanisms whose empirical value is, on present evidence, unproven at
matched wall time.

