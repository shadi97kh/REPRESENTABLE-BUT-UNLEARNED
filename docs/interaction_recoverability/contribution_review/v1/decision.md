# Contribution review decision

The bounded static review is complete. The combined mathematical result is ready for **external assessment**, with an unsent request and a self-contained evidence packet. It has not been externally reviewed. No backend, research, fitting or training cycle was opened.

The strongest supported claim is: in the specified five-law experiment, for known \(\theta\in[1/2,2]\), known \(0<\delta\le1/3{,}200{,}000\), the full original normalized analytic profile class and all four unknown coefficients in \([.5,1]\), the retained construction gives **pointwise regular attainment**, the worst-case efficient variance satisfies **\(\log(1+\sup V_{\rm eff})=\Theta(1/\delta)\)** uniformly in \(\theta\), and every balanced sample size \(N=5n\) obeys
\[
\mathcal R_{N,\delta,\theta}\ge\frac{e^{11/5}}{2{,}048{,}000}
\min\{1,e^{1/(640\delta)}/N,1/(\delta\sqrt N)\}.
\]
At \(N_\delta=5\lceil\delta^{-4}/5\rceil\), this implies risk at least \(c\delta\), whereas extending the historical exponential-profile upper order would imply \(O(\delta^2)\). Hence that rate cannot hold uniformly over the original full class with a \(\delta\)-independent constant. The fixed hard profiles could be disclosed to the estimator, so this does not isolate the cost of learning unknown profiles.

The [mathematical review](mathematical_review.md) checks the actual nonlinear likelihood path, global profile constraints, DQM domination, target derivative, allocation and testing factors. No fatal mathematical defect was identified in this static assessment. External scrutiny is still needed, especially of the constructive upper proof's growing empirical-quantile series, moving calibration queries, fixed charts, complete covariance, tail bounds and arithmetic realization. The lower proof does not rely on the upper proof; the paired efficient-variance order and pointwise attainment do.

The contribution is a concrete calibrated statistical example. Established adjoint-range, quantile, cohomological, finite-dimensional Riesz and two-point testing tools are not new principles. The exact overlap with DExtrI, including its broader appendix regimes, and the significance beyond this restrictive example remain open to external judgment. No publication-priority or acceptance claim is made.

| Question | Decision |
|---|---|
| Does a new learned GNN layer covered by this theorem exist? | **No.** The full-family quantile estimator is specified only. The executed neural-profile pilot has no graph; older RNA GNNs and the untrained graph CNP have different tasks and no theorem transfer. |
| Does the pilot justify a positive novelty claim? | **No.** Its saved 3.5569% aggregate reduction missed its frozen 10% screen; retain incremental-known and all negative findings. |
| What has siRNA demonstrated? | Data/provenance audits and ordinary predictor work, including an executed conditional RNA GNN that underperformed its guide baseline. For the combined theorem, siRNA remains a **motivating application only**. |
| Are calibrated private anchors, comparable replicated laws, compatible response model and held-out joint labels established? | **No.** No eligible six-condition panel was identified; exact missing evidence is listed in the scope report. |
| Is a finite-sample MSE upper bound, efficient estimator or true minimax rate proved? | **No.** A pointwise CLT and an influence-variance upper bound do not establish these. |
| Is exponential sample complexity for fixed accuracy proved? | **No.** The quadratic-path term prevents that inference from this lower bound. |
| Does the result cover unknown warps or historical \(\delta=.1\)? | **No.** The warp is known and the sufficient interval is explicit. |
| Historical reference efficiency? | **OPEN**, including the retained \(\delta=.1,\theta=1\) reference and efficient attainment. Worst-case small-\(\delta\) variance does not settle it. |
| Backend and practical route? | **NO-GO retained.** No useful finite-sample regime or fixed-resolution implementation is established, and no new implementation was authorized. |

Completed outputs:

- [Mathematical review and corollary](mathematical_review.md).
- [Actual architecture and siRNA evidence scope](gnn_and_sirna_scope.md).
- [External review request — draft, unsent](reviewer_request.md).
- [Dependency map](dependency_map.md), [manifest](manifest.json) and [packet instructions](README.md).
- `review_packet.zip` and its external `review_packet.zip.sha256` sidecar (distributed alongside this document).
- [Commands and preservation report](commands_and_preservation.md), [integrity checks](integrity_checks.json), and [preservation comparison](administrative/preservation_after.json).

All historical proofs, lower-bound v1/v2, frozen pilot inputs/predictions/configurations, negative results, ledgers and external_review/v2 were preserved within the documented hash scope. Only this fresh contribution_review/v1 directory was written. Administrative CPU/I/O had nonzero cost; measured helper-process costs and unmetered work are disclosed separately. No reviewer was contacted.
