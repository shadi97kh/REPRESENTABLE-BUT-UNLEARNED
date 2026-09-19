# Appendix completeness audit

The empirical continuation is `papers/interaction_recoverability_iclr2027/v2/`. Accepted v1 is preserved. The latest experimental brief expressly permits preserving mathematics in its separate draft; this empirical paper invokes no five-law theorem and does not conceal a theorem restriction in its appendix.

| Required dependency | Completed source |
|---|---|
| Complete notation, endpoint, known/missing inputs and sampling qualifications | `appendices/data.tex`, notation table |
| Raw acquisitions, pinned hashes, duplicate handling, exclusions and primary-table mappings | `data.tex`, generated source-hash/exclusion tables |
| Case-sensitive positional chemistry, linkage/conjugate/stereo conventions | `data.tex`, `model.tex`, full 45-name vocabulary |
| Grouping closure, exact seed and partition allocation, prior exposure | `data.tex`, generated source-membership table; final membership in supplement |
| Graph coordinates/edges/alignment, every feature and scaling constant | `model.tex`, context table and shift plot |
| Model forward map, parameter count, real gradient path | `model.tex`, main Equations 1–3, Proposition C.1 proof |
| Exact non-graph, no-chemistry, ridge/pairwise/CNN/tree controls | `model.tex`; invariance proof and numerical roundoff qualification |
| All losses, weights, optimizer settings, batches, stopping, seeds and selection | `training.tex`; all 20 development and 17 final fits |
| Checkpoints, resumption, target transform, prediction commitment | `training.tex`, Algorithm 1 dependencies |
| Every evaluation formula, tie convention, grouping/interval qualification | `evaluation.tex`; finite midrank-reference proof |
| All B1 aggregate/seed/study/pool results | `results.tex`, generated tables; 17 pools × 7 models |
| B2 eligibility, exact four-position chemical action and all 26 backgrounds | `data.tex`, `chemical_backgrounds.tex`, results/effect tables |
| B3 absence, SD-versus-SE and causal/theorem scope | `evaluation.tex`, `scope.tex` |
| Numerical precision, cost scaling, uncertified components | `model.tex`, `reproduction.tex` |
| Exact executed versus prepared commands, software, failures, cost boundaries | `reproduction.tex`, stage-attempt/software/fit tables |
| Prior primary-source comparisons and absence of invoked guarantees | `scope.tex`, citation audit |
| Every main figure/table/algorithm/equation dependency | `roadmap.tex` Table 5 (number verified in final build audit) |

No main-text claim terminates at an unproved internal theorem or a numerical surrogate for proof. Three elementary empirical/model propositions are proved in full. Bootstrap coverage is not claimed. Bulky exact rows, source copies and checkpoints complement the written specification; the appendix does not require reverse-engineering code to understand the method.

Missing scientific evidence is explicit: B3 four-condition labels, replicate design/covariance, independent study replication, complete training assay annotations, supported unseen-chemistry transfer, Davis strand mapping, and a new architectural contribution. These are not replaced by omitted proofs or appendix outlines.

Final PDF page boundaries, cross-reference checks, preserved-file hashes and rendered-page inspection are recorded in the final build audit.
