# Theorem-to-code map

| Mathematical step | Executed implementation | Evidence and limitation |
|---|---|---|
| Unique endpoint coordinates, linear observations and signed queries | `finite_range.py:exact_range` | Shared endpoint fixture; same endpoint coefficient cancels algebraically |
| Decoder equivalence with approximation error | Class-difference inequalities in `exact_range` | Constant-encoder, excluded-truth and eta fixtures; no learned eta guarantee |
| Compact polytope vertex enumeration | `solve_square`, active-set loop, rational feasibility checks | Complete nonlinear, missing-state and noisy exact solutions; proof in method_and_proof.md |
| Finite alternative encodings | `union_range` | Exact hull, empty/unresolved propagation fixtures; not adaptive membership certification |
| Floating-point/error handling | `rational`, explicit budget | Nonfinite rejection, no inferred zero for failures; only rational-problem certificate |
| Actual GNN forward map/no-message distinction | Preserved `sirna_gnn_empirical/v4/architecture.py:GraphModel.forward` | Executed 153,345-parameter and neighbor-substitution checks; not new architecture |
| Independent absent route | Actual architecture backward diagnostic | Zero absent-route data gradient, present-route nonzero gradient; scalar weight-decay counterexample |
| Both pair endpoint derivatives | Existing shared-weight loss; analytic torch diagnostic | +2 lambda v e and -2 lambda v e checked; no gradients through range inference |
| Actual training/selection | Preserved v4 `engine.py`, `train_focused.py`, frozen protocol and fit.json | 792 previous fits and 3,564 artifact hashes verified; zero new fits |
| Empirical metric/ensemble identities | `verify_reuse.py` | 36 score groups reproduced; seed-average decomposition checked |
| Promotion prerequisite | `runner.py:gate_block` | Pilot requires G1/G2 PASS; confirmation requires all three; missing runnable protocol also errors |

The exact finite-state result does not characterize all learned GNN responses, global neural optimization, or biological truth. G2 remains UNRESOLVED. No theorem is transferred from the unrelated five-law experiment.
