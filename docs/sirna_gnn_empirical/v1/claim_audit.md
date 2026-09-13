# Claim audit

| Claim | Evidence | Decision / scope |
|---|---|---|
| Actual learned chemistry GNN exists | `models.py`, profile, five main-GNN development/final fits and checkpoints; nonzero message gradients and changes | Supported. Conventional supervised message passing, 153,345 parameters. Two development fits + three final fits; profile is separate. |
| 4,766 admitted measured observations | `verified_eligibility.json`, row checks and observations | Supported: 2,927 ENsiRNA plus 1,839 primary-verified APP. No synthetic labels. |
| Valid frozen grouping | `split_manifest.jsonl`, `split_protocol.json` | Zero overlap of declared source/13-mer components. Does not establish biological independence or all forms of homology separation. |
| Predictions precede scoring | `prediction_commit.json`, `evaluation_results.json`, hashes | Supported; all 34,629 per-fit and 14,259 ensemble predictions are finite. No test-driven refits. |
| GNN beats validation-selected token CNN on APP MSE | .174726 vs .186286, paired difference -.011560, interval [-.018262,-.002634] | Supported conditionally within observed sources; not a new-study guarantee. |
| GNN beats strongest observed baseline | Tree .168597 vs GNN .174726; GNN-minus-tree interval [-.003419,.015526] | Not supported. Tree is strongest by observed test MSE, identified post hoc. |
| GNN produces useful APP candidate ranking | Top-five percentile .521442; within-assay rho -.003448 | No practical utility established. No prespecified practically useful gain threshold. |
| Chemical-effect transfer succeeds | 156 pairs, 26 exact backgrounds, 12 components; GNN pair MSE .292107 vs zero .280804 | B2 executed, effect-size improvement failed. All comparisons bundle four positions 6–9; no isolated GNA effect. |
| Measured two-change interaction evaluated | Maximum two chemistry states per matched block; zero complete rectangles | B3 incomplete. No four-condition labels; private anchors are not the obstacle for ordinary B3. |
| Biological mean-difference SE known | Source endpoint SDs, missing replicate n/shared-control covariance | Not supported. SD is not SE. |
| New GNN architecture or learning rule | Standard relation transforms/residual layers/readout; ENsiRNA and MEG-mod precedents | Not established. No novelty implied by implementation. |
| Five-law theorem covers GNN | Distinct inputs, supervision and estimator operations; absent anchors/latent-law verification | Not covered. Accepted theoretical manuscript remains separate, backend no-go and historical reference-efficiency question unchanged. |
| New wet-lab or causal validation | Public retrospective measurements only | Not performed; no causal mechanism, safety or clinical claims. |
| Numerical certification | Same-environment reload agreement, deterministic flags | Empirical reproducibility checks only, no rounding certificate/cross-hardware guarantee. |

Every numerical comparison is traced to the immutable current run. Historical pilot/reporter outputs are preserved and not reused as this campaign's evidence. Negative findings remain in all tables and plots.
