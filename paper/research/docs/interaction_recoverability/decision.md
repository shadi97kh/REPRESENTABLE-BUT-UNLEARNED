> Publication copy of `docs/interaction_recoverability/decision.md`. Mathematical and result content is retained; execution material is omitted. This dated record retains its original conclusions and limitations, including conclusions superseded by later stages. Original SHA-256: `26f9f33b77176be66d3919a484776b15bf38172e2e0267aa23d19d9b348ce4c8`.

# Interaction recoverability decision

**incremental-known.** The completed known-warp result has matching finite-sample rate orders for a selected interaction, but the estimator is ordinary two-quantile minimum-distance inversion. Its bilinear-completion control is the identical calculation. The broader unknown-mechanism contribution remains open, so no training authorization is requested.

**Strongest completed result — STANDARD RESULT APPLIED.** For the normalized class in [problem.md](problem.md), with known Gaussian noise, known exponential profiles/warps, known delta, fixed private anchors, four unknown bounded amplitudes and n independent outcomes at each of five observed actions,

\[
R_n(I_{12},\mathcal M_\delta)\asymp\min\{1,(n\delta^2)^{-1}\},\quad 0<\delta\le.1.
\]

Constants are uniform in n and delta across these separate known-delta classes, and are conservative. The [proof notes](proofs.md) derive the transformed-density score, bound its square integral by 65.6233873491624, propagate Hellinger affinity through the full product experiment, and give an amplitude-path lower bound matched in order by an empirical-quantile upper bound. No Wasserstein-to-testing inference is used. A secondary protected contrast equals a difference of observed means and has uniform O(1/n) MSE without recovering components. That is standard functional estimability, not a replacement primary endpoint.

**Exact remaining proof gap — OPEN.** There is no theorem here for learning unknown profiles/warps or separation from the same finite action dictionary while retaining a target-sensitive bound. An extension needs a specified nuisance class, a verifiable central-range/design condition and an estimator whose nuisance error is included, plus a distinction from existing inverse-functional/likelihood/completion theory. The restricted proofs contain no intentionally open lemma; they do not solve this broader step. An unrestricted eta*a1*a2 term is observationally invisible on all training axes and defeats robustness.

**Actual outputs — NUMERICALLY CHECKED.** The independent diagnostic reproduces the DOCX table and distinguishes joint W1 from the primary interaction gap:

| delta | Axis-2 W1 | Joint-action W1 | I12(T)-I12(Q) | Axis-2 Hellinger squared |
|---:|---:|---:|---:|---:|
| 0.1 | 0.556070894 | 3.295502812 | 2.776626586 | 0.000850339554 |
| 0.01 | 0.042849641 | 1.991464748 | 1.952692493 | 1.24270574e-05 |
| 0.001 | 0.004190336 | 1.896344042 | 1.892564740 | 1.29238818e-07 |
| 0.0001 | 0.000418112 | 1.887095558 | 1.886718583 | 1.29746526e-09 |

The limiting joint W1 and interaction gap are 1.886070833180237. Integration used [-12,12], absolute/relative requested tolerances 2e-11, explicit zero/crossing splits and transformed-density Jacobians. The largest W1 quadrature error estimate is 3.62e-11; the largest analytic omitted-tail upper bound is 2.21e-21. Hellinger inversion uses a 1e-12 root tolerance. Neither ordinary quadrature estimates nor floating-point evaluations of tail bounds are certified numerical enclosures. No samples, fitted estimators, bootstrap studies, tuning sweeps or benchmarks were run. The [diagnostic JSON](../../../../runs/interaction_recoverability/prep-20260912-082154/numerics-v2/diagnostics.json) retains every error estimate and analytic bound.

The [novelty matrix](novelty_matrix.md) checks the named primary results and direct adjacent collisions, including strongly identified functionals of weak nuisance functions. The DOCX and its nine embedded equation images were inspected. The separate `interaction_stability_check.py` was not found in the repository or attachment directory, despite follow-up; its original code remains unaudited. The independent harness reconstructs the equations supplied in both the prompt and DOCX. No original script was overwritten or falsely reported as executed.

The [estimator specification](estimator.md) gives the implemented, fitting-gated sample interface, deterministic algebra, costs and known/estimated/oracle information. The [siRNA requirements](sirna_requirements.md) document why current biological data do not validate this target. The [frozen unexecuted pilot](../../../../configs/interaction_recoverability_pilot.json) has SHA256 `1531f12364e099bc89be6bc7f06a8d7152e6af0f5a64437933bd89169031c8eb`. It includes equally informed minimum-distance/bilinear, ridge, likelihood and energy-score controls and labels native-method integration and oracle limits. The pilot is a falsification/control design, not an authorized run. Fitted runtime and optimizer behavior are unmeasured; an operation-count estimate is included instead of an invented complete compute cap.

**Resources.** The append-only ledger (local operational record; not published) belongs only to the new 600-CPU-second / 600-aggregate-job-wall-second allowance. The final snapshot (local operational record; not published) includes this report/check job after it exits, measured parent/waited-child CPU, conservative charge, separate task elapsed wall, memory and remaining allowance. Every metered job and descendant shared one enforced CPU affinity, with one numerical thread and no GPU. A remaining-budget watchdog kills the process group on timeout, with a clean-stop reserve. Initial text-only repository/prompt inspection preceded the meter and is explicitly unmeasured; 30 seconds are conservatively reserved/charged for it and launcher overhead. This is not presented as measured CPU. Unreaped descendants could be absent from measured counters, but inherited one-core elapsed time supplies the conservative bound for subsequent jobs. Remote retrieval was also charged where performed by the local wrapper, although the allowance permitted excluding it. Historical project balance remains unknown.

The verification record (local operational record; not published) confirms 521 original files unchanged, an empty tracked diff, CLI/help and syntax checks, and `14 passed in 0.31s`. The test count is implementation evidence, not the research outcome. The prior P4 and graph-CNP records remain unchanged.

**Next action:** subject the unknown-profile extension's proposed observation condition to mathematical falsification and the cited inverse-functional theory before spending a separate fitting budget. The current result does not justify a training request.
