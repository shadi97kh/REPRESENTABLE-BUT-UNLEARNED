# Workflow and forward-pass audit

The supplied architecture image is preserved through the verified v8 artwork and used unchanged on main-text page two. The final main figure caption is separate and below the image. No caption is embedded into a newly rendered artwork file. The diagram is a scientific schematic, not a new architecture or evidence that every accepted feature received a useful gradient.

| Artwork component | Actual implementation / qualification |
|---|---|
| 64 slots × 93 features | Frozen graph tensors, training-only active-node feature mask, ordered guide/passenger slots. |
| Linear 93→32 + ReLU | `GraphModel.embed`, followed by the padding mask. |
| Eight typed relations | Six directional backbone routes plus canonical/mismatch pairing. `A.sum((1,3))` supplies total receiver degree, floored at one. |
| Two message-passing layers | Per-relation neighbor aggregation, relation transform, biased self transform, ReLU, dropout .1, residual addition, LayerNorm, mask. |
| Shared support rule | Corrected backbone transforms reuse directional bases plus residuals gated by observed training support. Caption states this rule explicitly. No-message retains relation degrees and substitutes receiver for neighbor states. |
| Ordered H² readout | 64×32 states flatten into 2,048 coordinates; eight masked context coordinates give 2,056. |
| Dense 2,056→64→1 | ReLU and dropout .1 before the final linear normalized output. No sigmoid or clipping. |
| Activity scale restoration | Training equal-component mean and SD, with SD floor .05; no test recalibration. |
| Weighted activity supervision | Actual training partition component weights and standardized target. |
| Pair supervision | Separate B2 protocol; weighted squared standardized endpoint difference, shared model weights, selected loss multiplier. Both endpoint derivatives remain complete in the appendix. No new B2 training is needed for source robustness. |
| Fixed evaluation | Dropout disabled; endpoints and four B3 conditions evaluated with fixed checkpoint weights. |

The image's column-vector notation is the transpose of the appendix's row-vector equations. S7 and B3 visibility counts concern the actual frozen preprocessors, not the unmasked input schema. B3 exists in a separate primary panel; the APP inventory remains empty. The code still uses a conventional graph network with 153,345 scalars. No five-law theorem is transferred to it.
