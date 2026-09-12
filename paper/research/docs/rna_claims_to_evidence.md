> Publication copy of `docs/rna_claims_to_evidence.md`. Mathematical and result content is retained; execution material is omitted. This dated record retains its original conclusions and limitations, including conclusions superseded by later stages. Original SHA-256: `a74e2fec67f40e6a3ebb7d8bb973eb6a2c29b88b9cd9481ec94ca07902b55d24`.

# RNA application: claims to evidence

Generated 2026-09-05 by `scmp.rna.export_paper` from immutable result files. Task is **on-target siRNA efficacy and ranking**. No off-target or safety quantity is computed anywhere in `scmp/rna`.

## Supported

| claim | evidence | source |
|---|---|---|
| Huesken provenance preserved: records, guide length, unclipped labels, mapping coverage, historical split | 2431 records, all 21 nt, labels [0.000, 1.341] unclipped, 2431/2431 mapped, split [2182, 249] | data manifest |
| Assay identity taken from the primary source, not a redistribution | residual YFP fluorescence, H1299, 48 h, 34 constructs; 2 conflicts with the redistributed annotation recorded | data manifest |
| The historical split leaks; a sequence-cluster split does not | historical leak 87.1% by 13-mer cluster and 100% by gene; cluster split residual leak 0.00% | split audit |
| Ensemble machinery is exact where it claims to be | 102 checks, 0 failures: exact marginals match enumeration, every sample is a valid member, the deterministic mean bracket contains the exact mean, and yhat = b + a^T mu + gamma E[g] equals E_p[f] | ensemble audit |
| On-target efficacy is predictable from guide sequence alone under a leakage-free split | Spearman 0.6464 +/- 0.0160 over 5 seeds, P@20 0.130; sequence-cluster-disjoint split, single assay | evaluation |
| The u_i < l_j rule resolves a substantial fraction of pairs and certifies the MODEL ranking | 60/132 pairs resolved (45.5%); model ordering correct on 100% of them, as a certificate requires. Declared windows, NOT an RNA result | evaluation (engineering) |

## Not supported

| claim | status | reason |
|---|---|---|
| `mfe_gnn` arm | BLOCKED | missing: target_window, folding |
| `mean_adjacency_gnn` arm | BLOCKED | missing: target_window, marginals |
| `sampled_ensemble_gnn` arm | BLOCKED | missing: target_window, sampling |
| `anchor_A` arm | BLOCKED | missing: target_window |
| `fixed_A_plus_B` arm | BLOCKED | missing: target_window |
| `learned_A_plus_B` arm | BLOCKED | missing: target_window, policy |
| Structural sensitivity of efficacy (50/100/150 nt) | BLOCKED | the assayed context is a YFP reporter construct that is not in hand; native transcript windows would answer a different experiment |
| Learned-layer benefit | BLOCKED | requires the A+B arms |
| Useful residual nonlinearity on RNA | BLOCKED | requires the A+B arms; the ML-stage gate G4 already found the residual *hurts* held-out prediction (-0.268 Spearman) |
| Agreement of certified ranking with held-out measurement | BLOCKED | no measured ensemble readout exists for any window we can legitimately build |
| Chemical-variant prediction | BLOCKED | Huesken contains a single chemistry (unmodified); Davis 2025 has the chemistry field and is not yet acquired |
| Biological causation | NOT CLAIMED | the label is marginal over latent structures and is a reporter readout |

## Separation of claims

| kind | status |
|---|---|
| structural_sensitivity | distinct |
| chemical_variant_prediction | distinct |
| biological_causation | not_claimed |

`u_i < l_j` certifies that the **model** ranks i above j. It is not a statement about measurement. Agreement with held-out measurement is an empirical quantity, reported separately, and is currently not computable on any admissible RNA data.

## Gate status

| gate | status |
|---|---|
| ML G1 conditional vs unconditional | UNDETERMINED |
| ML G2 learned vs deterministic allocation | UNDETERMINED |
| ML G3 predictor noninferiority | PASS (baseline weak) |
| ML G4 residual predictive value | FAIL |
| DATA assay/construct identity | BLOCKING |
| DATA cross-dataset duplication | BLOCKING (datasets not acquired) |
| RNA pilot launch | REFUSED (no authorised budget) |
