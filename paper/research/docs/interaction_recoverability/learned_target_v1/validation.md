> Publication copy of `docs/interaction_recoverability/learned_target_v1/validation.md`. Mathematical and result content is retained; execution material is omitted. This dated record retains its original conclusions and limitations, including conclusions superseded by later stages. Original SHA-256: `38690b4ba9ed7b2560d8ddd9929d94e04b56b1318757c9076f060f7563ab5f83`.

# Correctness and execution validation

Status: eight focused correctness tests passed before the frozen pilot. See the machine-readable [validation output](../../../../results/runs/interaction_recoverability/learned-target-v1-20260912-164000/validation.json). This count is for the new learned-target namespace, not the historical suites.

| Consequential risk | Completed check |
|---|---|
| Wrong likelihood gradient at fixed outcome | Independent transformation-score derivative and central finite differences, including all shifts and both profile blocks |
| Non-normalized density or nonzero-mean score | Independent log-outcome density integration and score integration, with Gaussian-tail treatment |
| Invalid profile/perturbation constraints | Projection and normalization fixtures; global analytic envelope proofs in `method.md`; finite admissible path lengths |
| Missing target-gradient terms | All four coefficient and both profile contributions compared to independent finite differences |
| Incorrect inverse/ridge algebra | Linear Gaussian inverse fixture with exact Gram, target bias and variance |
| Suppressed unidentified target direction | Score-null target fixture retains residual bias; well-conditioned comparison and coupled zero-delta cancellation diagnostic |
| Incorrect nonlinear bias accounting | Finite-path score identity and a nonzero nonlinear remainder; bias sign and variance decomposition |
| Incorrect allocation or leakage | Deliberately unequal stratum weights/counts; disjoint fold indices; strict source archive keys; no evaluator import in learner |

The recorded CLI validation completed in 1.92 seconds with eight passes and one PyTorch warning about converting a requires-grad tensor to a scalar in a diagnostic. The earlier direct focused pytest invocation also passed eight tests. The warning does not skip the derivative checks; it concerns scalar logging. It is retained in the captured output.

The development execution covered both fitted profile families used as **estimators**, all three F/C/E rotations, all ten controls, endpoint searches, saved predictions and separate development truth evaluation. That is one development dataset, not another final replicate. Timing selected the predeclared full plan before final sample generation. Source loss, iteration/checkpoint choice, integration differences, model payloads and endpoint search traces are saved.

The live budget watchdog was tested with a one-second job spawning a sleeping child in a new process group. The wrapper returned the expected timeout status 124; the child was subsequently absent from `/proc`. The intentional timeout and its partial PID artifact remain recorded. This test demonstrates the observed cleanup path, not protection against malicious unobservable process/affinity escapes. Admitted jobs inherit one-core affinity and one numerical thread; the watchdog checks process-tree affinity and terminates observed descendants.

The weighted profile-norm quadrature has a 4096-versus-2048 spectral-norm discrepancy of **1.368362873580413**; smallest M eigenvalue is **0.078445254127133**. It is not a certified integral enclosure. Development fold-zero score-Gram refinement differences were `1.2897e-7` (neural) and `1.1900e-7` (basis), target-gradient differences `7.8164e-10` and `9.4118e-10`, and centering maxima `9.5722e-9` and `8.4524e-9`. These are actual development diagnostics, not final maxima or error guarantees. Final numerical diagnostics are aggregated in `results.md`.

No historical benchmark, P4 training, old preparation suite or instability grid was rerun. Preservation checks compare hashes and git state. Final scoring checks source hashes and freezes all prediction hashes before opening truth. Tests support code correctness on the checked cases; they do not prove a useful nonlinear risk rate, optimization optimality, finite-sieve completeness, biological assumptions, or a statistically reliable pilot win.
