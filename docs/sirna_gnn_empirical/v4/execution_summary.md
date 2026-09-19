# Focused robustness results — manuscript v9 / empirical v4

| Sensitivity        | Method         |   Grouped MSE |   APP MSE |   S7 MSE |
|:-------------------|:---------------|--------------:|----------:|---------:|
| source_excluded    | GNN            |      0.064083 |  0.163060 | 0.119091 |
| source_excluded    | No-message     |      0.062518 |  0.163657 | 0.122147 |
| source_excluded    | Chemistry tree |      0.052931 |  0.165515 | 0.097356 |
| source_excluded    | Token CNN      |      0.062798 |  0.162649 | 0.110302 |
| primary_reanchored | GNN            |      0.074727 |  0.162821 | 0.105208 |
| primary_reanchored | No-message     |      0.074554 |  0.162397 | 0.107330 |
| primary_reanchored | Chemistry tree |      0.069160 |  0.189936 | 0.093435 |
| primary_reanchored | Token CNN      |      0.075276 |  0.164069 | 0.101890 |

Executed **792 new fits**: 432 development and 360 final; **0 failed fits**. All final stochastic blocks retain ten seeds. Executed neural updates: **476,867**; selected-checkpoint updates: **166,067**. Forty original final checkpoints were verified and reused; B2/B3 original fits and predictions remain separate, without new fitting.

- source_excluded, ENsiRNA_grouped: GNN minus matched no-message MSE **+0.001565 [-0.000851, +0.004494]**, conditional on fixed ensembles and observed components.
- source_excluded, APP: GNN minus matched no-message MSE **-0.000597 [-0.001670, +0.000338]**, conditional on fixed ensembles and observed components.
- source_excluded, Davis_S7: GNN minus matched no-message MSE **-0.003057 [-0.004147, -0.001774]**, conditional on fixed ensembles and observed components.
- primary_reanchored, ENsiRNA_grouped: GNN minus matched no-message MSE **+0.000173 [-0.000838, +0.001290]**, conditional on fixed ensembles and observed components.
- primary_reanchored, APP: GNN minus matched no-message MSE **+0.000424 [-0.000824, +0.001550]**, conditional on fixed ensembles and observed components.
- primary_reanchored, Davis_S7: GNN minus matched no-message MSE **-0.002122 [-0.003151, -0.000894]**, conditional on fixed ensembles and observed components.

Both new S7 sensitivities establish a narrow conditional ensemble improvement over matched no-message; tree and CNN nevertheless have lower observed MSE. A general GNN advantage is not established by these heterogeneous sensitivities. Original APP tree–GNN differences are dominated by squared mean bias (99.53%). The original S7 routing improvement is a narrow conditional positive result; it does not establish a message-passing advantage. B2 sequence-specific measured chemistry prediction remains unestablished, and original B3 does not beat zero interaction. These are retrospective measured-data analyses, not independent external or new wet-lab validation.

Primary-source adjudication reproduces 299 residuals above .05 among 1,604 linked reported identities (max 1.0852170325). The 1,924 admitted Bramsen rows partition into 1,604 linked, 201 unmatched, 76 duplicate release mappings and 43 without a unique primary measurement. The separate primary-assay sensitivity contains 2,607 rows and preserves all released values; it is not a certified release repair. Full molecular/assay traceability remains incomplete. The primary article reports triplicate assays repeated twice, but shared-control covariance is unavailable. B2 replicate counts remain unresolved.

Exact input counts, chemistry-name support, all 140 cancellation flags, dependent leave-margin and reference-shift checks are in `visibility/`. Cross-campaign ranking uses identical 1,839 records and 17 pools; corrected exact tie handling retains the principal reversal. Earlier constants and some chemistry-invariant scores change under the common definition, and their originals are preserved. Published ENsiRNA/MEG exact reproductions remain blocked (zero fits); acquired Davis S1 remains quarantined.

Main pages **1–9**; statements **10**; references **11–14**; appendix **15–53 (39 pages)**. The workflow is on page **2**. Source validation, external build, independent ZIP build, proof preservation and all-page rendering are recorded separately.

Deliverables:

- PDF: `papers/interaction_recoverability_iclr2027/v9/artifacts/main.pdf`
- Four-entry source ZIP: `papers/interaction_recoverability_iclr2027/v9/artifacts/manuscript_source.zip`
- External official dependencies: `papers/interaction_recoverability_iclr2027/v9/artifacts/typesetting_dependencies.zip`
- Final referenced figures: `papers/interaction_recoverability_iclr2027/v9/source/figures/`
- Full source adjudication: `runs/sirna_gnn_empirical/v4-20260914T071926Z/adjudication/row_reconciliation.csv`
- Predictions and checkpoint directories: `runs/sirna_gnn_empirical/v4-20260914T071926Z/evaluation`, `runs/sirna_gnn_empirical/v4-20260914T071926Z/fits`
- Separate indexed evidence archive: `runs/sirna_gnn_empirical/v4-20260914T071926Z/delivery/robustness_code_results.zip`

Recorded completed stage/diagnostic costs so far: 4661.377 child CPU seconds and 4448.602 aggregate command wall seconds; max reported child RSS 3228.1 MiB. GPU fitting-region elapsed 3501.493s overlaps stage time and is not GPU kernel time or energy. Direct administrative work has additional unmetered nonzero cost; build/package receipts are separate. Historical datasets, proofs, predictions, configurations, checkpoints, logs and negative results remain preserved. No push, publishing or contact occurred.

The named handoff ZIP and review files were absent from the accessible workspace and attachments; the requested local path has not been supplied. The concrete review corrections in the root instruction were independently checked. This delivery does not claim extraction or inspection of missing files.
