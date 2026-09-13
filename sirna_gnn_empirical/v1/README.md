# Executed measured-siRNA benchmark

This is an actual trained conventional chemistry-aware nucleotide GNN and six controls, evaluated against measured siRNA outcomes. It is not the five-law quantile estimator and makes no architecture novelty claim.

The active immutable campaign path is in `active_run.txt`. Run `bash sirna_gnn_empirical/v1/run_all.sh` to verify completed stage hashes and resume any unfinished stage. For an independent replay, use `--run runs/sirna_gnn_empirical/NEW-UNIQUE-ID`; do not overwrite existing evidence. `--stage STAGE` and `--through STAGE` are supported. The full source evidence JSON is required at repository root; dependencies and public source availability are specified in the manuscript appendix.

Stages: acquire, dependencies, inspect_sources, supplements, patent_sources, primary_tables, reconcile, verify_assays, split, graphs, profile, develop, freeze, train, predict, evaluate, figures, paper. Profile uses measured training labels. The plan freezes after timing. Fits save model/optimizer states after each epoch. Predictions are hashed before test labels are scored. No post-test tuning is supported.

Final measured population: 4,766 rows; train 2,626, validation 103, internal test 198, APP test 1,839. Seven methods, 17 final fits, 6,636 final neural updates (1,701 main GNN), plus preserved development/profile work. All final stochastic models use seeds 1103/2207/3301.

APP GNN MSE .174726 versus chemistry tree .168597: no improvement over the strongest baseline is established. B2 ran on 156 four-position chemistry-pattern comparisons across 26 backgrounds/12 sequence components; GNN effect MSE .292107 versus zero-difference control .280804. B3 has zero eligible four-condition sets and remains incomplete. Results are retrospective measurements, not new wet-lab validation.

Complete claim/citation/appendix audits and execution summary are in `docs/sirna_gnn_empirical/v1/`. The empirical paper continues in `papers/interaction_recoverability_iclr2027/v2/`, preserving v1 and historical proofs. Original machine logs are retained locally; the complementary scientific ZIP omits machine-specific logs for anonymity. Do not interpret optimization seed SD or conditional cluster intervals as biological replicate uncertainty.
