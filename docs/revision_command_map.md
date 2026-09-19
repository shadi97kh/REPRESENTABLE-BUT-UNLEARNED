# Revision command map

Every `scmp.revision.*` interface requested by the action plan, with its actual
implementation status. **Requested** means the plan asked for it; **implemented**
means it exists here and has been run; **planned** means it is not written yet
and must not be invoked.

All commands support `--help`. Every compute command takes `--dry-run` and prints
a budget report before doing work. Nothing under `runs/` that a previous run
produced is ever modified: corrections are written as new records that cite the
record they correct.

## Implemented and run

| Requested in plan | Implemented as | Status | Writes |
|---|---|---|---|
| — (master instruction: record source hashes) | `python -m scmp.revision.provenance` | implemented | `runs/revision/snapshot-*.json` |
| `scmp.revision.audit_report --report run_report.pdf --runs runs` | `python -m scmp.revision.audit_report --report docs/run_report.md --runs runs` | implemented | `runs/revision/audit_report-*.json` |
| `scmp.revision.audit_math --exact-arithmetic` | `python -m scmp.revision.audit_math --exact-arithmetic` | implemented | `runs/revision/audit_math-*.json` |
| — (plan: "a machine-readable corrected ledger") | `python -m scmp.revision.ledger` | implemented | `docs/corrected_ledger.json` |
| — (P0 traceability rescue) | `python -m scmp.revision.rerun_unrecorded [--dry-run] [--reparse REC]` | implemented | `runs/revision/rerun_unrecorded-*.json` |
| `scmp.revision.predictor_pilot --config configs/revision_predictor.yaml [--dry-run]` | same, plus `--stages S0..S5` | implemented | `runs/revision/predictor_pilot-*.json` |
| `scmp.revision.profile_verifier --config configs/revision_profile.yaml` | same, plus `--dry-run` | implemented | `runs/revision/profile_verifier-*.json` |
| `scmp.revision.verifier_pilot --config configs/revision_verifier.yaml [--dry-run]` | same | implemented | `runs/revision/verifier_pilot-*.json` |
| `scmp.revision.rna_audit --config configs/revision_rna_data.yaml` | same, plus `--skip-insert` | implemented | `runs/revision/rna_audit-*.json` |
| `scmp.revision.rna_pilot --config configs/revision_rna.yaml [--dry-run]` | same, plus `--limit` | implemented | `runs/revision/rna_pilot-*.json` |

### Deviations from the requested signature

- `audit_report` takes `--report docs/run_report.md`, not the PDF. The PDF is a
  rendering of that Markdown (`docs/run_report.pdf`, same content, generated from
  it). The path is **recorded, not parsed**: no number in the audit is scraped
  from the report. Every figure is recomputed from the immutable record, which is
  what "trace to immutable records" requires. Parsing the prose would make the
  report its own evidence.
- `audit_math` runs its Fraction-based checks by default; `--exact-arithmetic` is
  accepted and is a no-op flag retained for the plan's exact invocation.

## Planned — not implemented, do not invoke

| Requested in plan | Phase | Blocking dependency |
|---|---|---|
| `scmp.revision.audit_numerics --config configs/revision_numerics.yaml` | P4 | operator scope fixed |
| `scmp.revision.check_certificates --config configs/revision_enclosed.yaml` | P4 | `audit_numerics` |
| `scmp.revision.freeze_protocol --config configs/protocol_v2.yaml` | P6 | P0–P5 |
| `scmp.revision.evaluate --config configs/protocol_v2.yaml [--dry-run]` | P6 | frozen protocol v2 |
| `scmp.revision.export --config configs/protocol_v2.yaml` | P6 | evaluation records |

## Existing commands reused unchanged

These already exist and are correct; the revision cycle calls them rather than
reimplementing them. Preserving their filenames is deliberate.

| Command | Purpose |
|---|---|
| `python -m scmp.audit_bounds` | exhaustive bound-soundness audit |
| `python -m scmp.audit_conditional` | conditional-support soundness |
| `python -m scmp.audit_oracles` | oracle correctness vs independent enumerator |
| `python -m scmp.check_certificates` | independent certificate checker |
| `python -m scmp.evaluate_gates` | protocol v1 gates (historical; not re-scored) |
| `python -m scmp.reproduce` | hash-checked reproduction |
| `python -m scmp.export_results` | tables from immutable hashed files |
| `python -m scmp.rna.audit_data` / `audit_splits` / `audit_ensemble` | RNA audits |

## Test entry points

```
python -m pytest tests/revision/ -q      # P0 regression tests (24)
python -m pytest -q                      # full suite
```
