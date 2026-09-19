# Appendix completeness audit

| Main dependency | Complete material | Verification coverage |
|---|---|---|
| Main Figure 1 / architecture | appendices/architecture.tex; appendices/training.tex | Forward tensors, explicit initialization and count, masks, fallback, no-message proof, exact normalized losses and both endpoint derivatives. |
| Main Table 1 and support figure | appendices/data.tex; complete support tables | Source/admission/grouping and all partition relation counts, rare/absent cases, per-layer gradients/deltas and reachable-scalar upper bounds. |
| Activity table and figure | appendices/metrics.tex; appendices/results.tex | All cohort/arm/weighting/seed/fold/assay error and calibration panels; trained/oracle constants. |
| Direct comparison table | appendices/metrics.tex; complete direct/seed tables | Matched differences, all 10,000-resample targets, training variation separate; exact ensemble decomposition. |
| B2/ranking table and figure | appendices/data.tex; training.tex; metrics.tex; results.tex | Every pair identity/SD, objective/folds, control definitions, endpoint accuracy, exact ties and all pool overlap. |
| B3 table | appendices/published_methods.tex; complete B3 tables | All states, aliases, four measured values and SDs, source holdout, common-reference limitations. |
| Main optimization table and all executed fits/selection | appendices/training.tex; all-fit and selection tables | Every fit ID has a closed exact path/config dictionary; all candidate selections and outcomes, not just winners. |
| Reproduction/resources/failures | appendices/reproduction.tex | Arithmetic, versions, exact command snapshot, resume qualifications and actual cost scopes. |
| Retained mathematical claims | appendices/architecture.tex; appendices/metrics.tex | Full finite-model proofs; no appeal to an omitted internal theorem. |
| Prior work | appendices/published_methods.tex; citation_audit.md | Primary source/claim checks; exact attempted published reproduction and blockers; no imported published GNN risk theorem. |

The unrestricted appendix contains current evidence and model-specific derivations, not reprints of superseded papers. It does not end a required proof at a repository pointer. CSVs preserve full numerical precision in addition to the printed synchronized panels. Separate build and visual audits record page boundaries and complete rendered-page inspection.
