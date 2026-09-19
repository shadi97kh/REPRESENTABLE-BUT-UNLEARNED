# Reporting erratum: efficiency solve-residual scope

This is **new commentary**. The historical efficiency argument, code and raw results remain unchanged. Verification used only the saved JSON values and document text; no solve, integration, test or diagnostic was rerun.

## Original wording

Original repository file:
`docs/interaction_recoverability/efficiency_audit/v1/efficiency_argument.md`, line 226, in “Deterministic values and numerical qualifications”:

> Maximum solve residual is \(7.24\,10^{-15}\); discarded-\(b\) norm is zero.

The paragraph follows the finer-configuration rank table, but does not explicitly qualify that maximum by configuration. Read as a maximum across both saved configurations, the wording is incorrect.

## Correct scope and values

**Corrected statement:** The maximum solve residual in the **fine configuration** is
\(7.23776861244768\times10^{-15}\), approximately \(7.24\times10^{-15}\).
Across **both saved configurations**, the maximum is
\(9.540492121253299\times10^{-15}\), approximately \(9.54\times10^{-15}\), in the coarse twelve-direction record.

Original raw source:
`docs/interaction_recoverability/efficiency_audit/v1/raw.json`.

| Configuration | Directions | Exact saved solve_residual | JSON pointer |
|---|---:|---:|---|
| coarse | 4 | 2.220446049250313e-15 | `/results/0/score_geometry/lower_bounds/0/solve_residual` |
| coarse | 8 | 5.5502823933978995e-15 | `/results/0/score_geometry/lower_bounds/1/solve_residual` |
| coarse | 12 | **9.540492121253299e-15** | `/results/0/score_geometry/lower_bounds/2/solve_residual` |
| fine | 4 | 2.220446049250313e-15 | `/results/1/score_geometry/lower_bounds/0/solve_residual` |
| fine | 8 | **7.23776861244768e-15** | `/results/1/score_geometry/lower_bounds/1/solve_residual` |
| fine | 12 | 6.202610464077155e-15 | `/results/1/score_geometry/lower_bounds/2/solve_residual` |

These are literal saved-field extracts, not newly calculated residuals. The manifest records the original raw-file path and SHA256; the full Gram arrays are deliberately omitted from this reading packet.

## Consequence

This corrects the scope of a reporting maximum. It changes no saved matrix, coefficient, target derivative, covariance, submodel lower-bound value or corrected influence value. Both configurations retained all prescribed eigenmodes. Small linear-solve residuals still do not certify integration accuracy, efficient variance or a feasible efficient estimator.

The approximate 38,776-fold efficiency gap remains unresolved. The unfavorable influence scale, 0.081285% median reduction and backend no-go are unchanged. No performance or novelty claim is strengthened.

See the historical [efficiency argument](evidence/efficiency_audit/v1/efficiency_argument.md), the new [review brief](review_brief.md), and [claim map](claim_map.md). The nearby pseudoinverse-agreement statement also sits with the fine-configuration table; the package does not reinterpret it as an all-configuration certificate.

