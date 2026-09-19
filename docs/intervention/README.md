# Intervention-explainer preparation — September 12, 2026

**Decision: incremental known method.** The requested preparation is complete. The proposed response-only architecture is a graph conditional neural process, and the generic CNP comparator is the identical class. Pair and matched-output losses are implemented as ablations; their combination does not establish novelty. No predictor or explainer training ran, and there is no positive held-out evidence.

The revised question is whether those controlled losses improve adaptation to new predictor functions over the same response-only operator. That is an empirical question about an incremental method. A substantive-contribution route needs a concrete mechanism beyond the identified equivalence before training is justified. The [source-backed matrix](novelty_matrix.md) records the remaining differences, code availability and bounded search through September 12, including the unresolved identity of “grammar-guided local data.”

The implementation passes **18 correctness checks**. It reconstructs graph/descriptor/stereo inputs, preserves identity centering and query-set invariance, conditions on observed responses, exposes the three loss arms and strong direct controls, and exports an [80-row graph/model/action response table](../../runs/intervention/preparation-v2/responses.csv). Those rows are fixed analytic model-query fixtures, not biological measurements or held-out results. Exact evaluation of all eight states is already faster than the untrained operator's forward pass in the measured tiny fixture; no efficiency break-even exists in that observed regime. The analytic Gaussian function family also admits a correctly specified GP conditional-mean control.

The data audit reconciles 1,200 AGILE source-workbook values and produces [3,500 provenance rows](../../runs/intervention/audit-v2/rows.jsonl). AGILE/LANTERN same-code labels agree, while 100 product representations differ stereochemically; chemical joins alone misleadingly produce 100 label conflicts. TransMA has conflicting repeated structures within its splits. LANTERN's pinned `Murcko_scaffold` file shares six train/test scaffolds, whereas `scaffold_balanced` shares none. Mixture identity, per-row replicates and checkpoint training membership remain limits on a biological study.

| Reviewable output | Location |
|---|---|
| Novelty, exact estimands/access/inputs/targets/costs and source limits | [Matrix](novelty_matrix.md), [structured version](novelty_matrix.json) |
| Mathematical definition, transformations, operator, controls and information-flow checks | [Method specification](method.md), [implementation](../../intervention/operator.py) |
| Pinned sources, experimental provenance, stereo/conflict/split/checkpoint results | [Data report](data_audit.md), [source manifest](../../data/intervention_sources/manifest-v1.json), [supplement manifest](../../data/intervention_sources/supplement/manifest-v1.json) |
| Frozen primary endpoint, splits, baselines, cost allocation and stop rules | [Configuration v2](../../configs/intervention_screen_v2.json), [SHA256](../../configs/intervention_screen_v2.json.sha256) |
| Actual diagnostic costs, partial prospective estimate and missing measurements | [Resources](resources.md), [machine-readable summary](../../runs/intervention/resource_summary-v1.json) |
| Final checks and preservation evidence | [18 passing tests](../../runs/intervention/tests-v4.json), [integrity record](../../runs/intervention/integrity-v1.json), [bundle hashes](../../runs/intervention/bundle-v1.json) |

Configuration v2 supersedes v1 before any training or held-out execution. It adds the known-prior GP/tree controls and specifies exact source RNG, empty-query handling and initialization reduction; the endpoint, time interval and continuation threshold are unchanged. Both frozen files and earlier audit/test/profile versions remain available. No external PGExplainer/FastSHAP/GraphSHAP-IQ/SPEX reproduction, published delivery-model inference, RNA folding pipeline or end-to-end trainer is claimed. Applicability and missing integration are explicit in the method and configuration.

The standalone checks can be reproduced from the repository root:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python -m pytest -q tests/intervention
```

The completed final profile was generated with the following preparation-only command. Its existing destination is exclusive: a deliberate repeat needs a fresh output directory.

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python -m intervention.prepare --out runs/intervention/preparation-v2
```

The completed negative [Grammar-CertMP/P4 decision](../p4_actual_model_executed-20260911-171200-784773-final.md) and its recorded artifact hashes are preserved. Its training route was not reopened. The old remaining resource allowance is still unknown. The current authorization covered preparation and small checks; no new training allowance has been inferred or requested for a mechanism already equivalent to prior work.
