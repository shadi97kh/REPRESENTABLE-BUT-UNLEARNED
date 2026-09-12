> Publication copy of `docs/interaction_recoverability/sirna_requirements.md`. Mathematical and result content is retained; execution material is omitted. This dated record retains its original conclusions and limitations, including conclusions superseded by later stages. Original SHA-256: `64161fc27b87f1b3243bd2c3e87d2dcb34befd785f2b7e0846e36f8ca5872d40`.

# siRNA evidence requirements and scope

The mathematical result concerns independent outcome distributions at declared actions on a fixed background. It does not validate siRNA efficacy, modification interactions, a causal molecular mechanism or explanations of a fitted GNN. A GNN is unnecessary for the completed known-basis theorem.

| Required field/design property | Why it matters | Current evidence limit |
|---|---|---|
| Exact guide and passenger sequences, target background and strand orientation | Defines the unchanged background and position correspondence | Existing reporter work includes unresolved construct/isoform context; different sequences are not replicates |
| Exact chemical modification, stereochemistry, linkage and position for each action | Defines a feasible categorical edit and its identity | An unmodified sequence table cannot support unseen-chemistry claims |
| Reference, both single edits and specified joint edit; additional anchors if assumed | Defines I12 and permits held-out evaluation of missing-combination transfer | Obtaining a supplementary file does not establish this coverage |
| Matching assay, delivery, cell system, dose and time | Makes observed distributions and means comparable | Pooling assay groups would change the observation law |
| Declared endpoint scale and normalization, including raw controls | Interaction depends on scale; shared controls induce dependence | Log means, mean logs and ratios of means are different targets |
| Independent biological replicates per action, batch/well IDs and technical-repeat nesting | Justifies the stratified product law and quantile learning | One efficacy value per sequence does not give a replicated response distribution |
| Raw values, censoring/missingness, source workbook/row, version and hash | Preserves provenance and limits selection bias | Missing or selected observations require an explicit sampling model |
| Split units grouping shared background, molecules and controls | Prevents pair overlap from inflating evidence | Many rectangles from one background are not independent test cases |
| Calibration of any known warp, private-anchor effect or latent-noise assumption | Separates a usable method from an oracle calculation | None of these structural calibrations is established for current siRNA data |

The existing [P3 audit](../p3_rna_audit_and_pilot.md) describes reporter efficacy, unresolved or inferred construct context and missing modified-chemistry evidence. Its prior negative result is preserved. The [P4 actual-model decision](../p4_actual_model_executed-20260911-171200-784773-final.md) also remains negative. The [delivery-lipid audit](../intervention/data_audit.md) concerns an mRNA delivery material with mixture, replicate and representation limitations; it is not an siRNA modification factorial experiment. These summaries were read, without new broad biological acquisition or processing.

A suitable later benchmark would have a fixed sequence background and two named modifications, independent batch-resolved measurements at the reference and each single edit, and independently measured joint outcomes withheld from estimation. Multiple backgrounds would create separate independent evaluation units with prespecified grouping. Distributional methods need enough replicates within each actual stratum, not artificial pooling across backgrounds. A theorem requiring a continuously adjustable action and interval coverage would need real dose/strength sweeps; binary modification indicators do not substitute for them.

For a measured endpoint, the explanation is an estimated feasible-change contrast with its sampling and model assumptions. For a frozen predictor, the explanation is the difference between two model queries after rebuilding every dependent input. Predictive fidelity alone does not imply agreement with measured biological effects. GraphSHAP-IQ, SPEX, direct querying and a graph encoder remain known components; none repairs absent measurements or establishes the new statistical contribution.

New combinations of already observed modifications, new sequence backgrounds, new assay contexts and unseen modification chemistry are separate generalization axes. The present work supplies synthetic mathematics and deterministic fixtures only. Biological validation is currently unsupported, and the proposed fitted control pilot therefore contains no siRNA efficacy claim.
