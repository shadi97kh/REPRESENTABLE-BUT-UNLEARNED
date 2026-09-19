# Workspace cleanup — 14 September 2026

Removed 3,210 unnecessary generated files, totaling 751.3 MiB of file contents. Actual filesystem allocation recovered may differ.

Removed: temporary full-page PDF extracts, raster inspection previews/contact sheets, byte-identical source-ZIP extraction copies, generated independent-build artifacts, project Python bytecode/pytest caches, and untracked TeX intermediates. Exact removed paths, previous hashes and reasons are in `deletion_manifest.json`.

Preserved: all tracked file bytes and existing tracked modifications; every version's final manuscript PDF and source ZIP; mathematical source/proof files; scientific code, datasets, results, checkpoints, ledgers and logs; review archives and evidence packages; task briefs and attachments, including the files open in the IDE; and dependencies still used by the scientific code or renderer. No content under `runs/`, `data/`, `review_packages/` or `docs/interaction_recoverability/` was selected for deletion.

Validation: 1,750 tracked files and 634 protected source/deliverable files retained their hashes. The scientific publication verifier passed for 1,330 current local artifacts. All seven manuscript-v8 delivery artifacts match their delivery-manifest hashes. No new fits, inference, scientific analyses or historical ledger writes occurred. Nothing was pushed for this cleanup because no tracked file changed.

Earlier preservation and rendering manifests remain historical records of their original checks. They may name disposable intermediate files that this cleanup now removes; this dated deletion manifest accounts for those removals. Source-build success and visual-audit reports are retained. Regenerate temporary output from the retained source/PDF when needed.

Administrative copying/inspection/hashing and deletion incurred nonzero cost. The cleanup script's measured wall time was 6.527 seconds; surrounding inspection and verification time was additional and unmetered.
