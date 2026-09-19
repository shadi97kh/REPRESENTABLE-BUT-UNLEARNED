# Formatting changes: v2

The v1 directory and its ZIP are unchanged.

## One commentary correction

In claim_map.md, the Bennett source-condition sentence contained the unprotected expressions Range(P*P) and Range(P*). The installed Markdown renderer treated the two asterisks as an emphasis pair spanning the intervening words.

Both existing expressions are now enclosed in inline-code delimiters: `Range(P*P)` and `Range(P*)`. Their characters and mathematical meaning are unchanged. No global replacement was used.

The [review brief](review_brief.md) and [erratum](errata.md) are byte-for-byte unchanged. Their expansion, observation/influence indices, contraction statement, precision wording and numerical qualifications required no correction. In particular, 7.24e-15 remains the fine-configuration maximum and 9.54e-15 the maximum across both configurations.

All 25 historical evidence copies are byte-for-byte identical to v1, including its explicit labels for omitted repository references. No theorem, assumption, constant, quantifier, estimator definition, influence identity or scientific conclusion changed.

## Packaging and verification

The manifest records this commentary change and its new hash. The ZIP was rebuilt for the new version and extracted into a temporary directory for active-link and hash inspection. A separate review_packet.zip.sha256 sidecar records the new archive hash and is excluded from the archive.

The existing markdown-it renderer confirmed the accidental emphasis and the corrected inline-code rendering. Markdown table structure was inspected. Mathematical expressions were checked in source text against the retained statements; no math-capable visual preview was available. Plain Markdown rendering is not claimed as visual verification of LaTeX.

Only document inspection, copying, rendering, compression and integrity hashing were performed. These have nonzero CPU cost. No scientific test, numerical calculation, integration, sample, fit, training, simulation or benchmark was run. Historical evidence, v1 and every compute ledger were preserved.

Scientific status remains: local attainment retained; efficiency unresolved; backend no-go retained; publication novelty unestablished. This remains a mathematical reading packet, not a runnable reproduction environment.
