# Contribution review packet

Start with [decision.md](decision.md), [mathematical_review.md](mathematical_review.md) and [gnn_and_sirna_scope.md](gnn_and_sirna_scope.md). The [review request](reviewer_request.md) is a draft and has not been sent. This is internal review preparation, not external mathematical certification.

`evidence/` retains original repository-relative paths. Every included evidence file is an unchanged byte copy. The [dependency map](dependency_map.md) explains the necessary mathematical chain, historical comparisons, source-code/run evidence and omitted execution dependencies. Read copied code as text; it is not an executable reproduction package, and no scientific program is needed to assess the review.

[manifest.json](manifest.json) records SHA-256 and byte length for payload files. It excludes itself and final packaging reports to avoid self-reference. [integrity_checks.json](integrity_checks.json) checks paths, bytes and local links, with its exact scope stated. The outside `review_packet.zip.sha256` sidecar authenticates `review_packet.zip`; neither file is inside the archive itself. All proof, review and evidence navigation is checked inside the archive.

[commands_and_preservation.md](commands_and_preservation.md) reports permitted administrative work, measured helper costs, unmetered costs, and preservation checks. Existing proofs, predictions, configuration, ledgers and the earlier external_review/v2 packet were not edited.
