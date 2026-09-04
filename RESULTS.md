# Results -- first full run

Run date 2026-09-04. Pre-registration frozen at tag `prereg` before any experiment ran.
Environment: numpy 2.4.2, ViennaRNA 2.6.4, torch 2.11.0. Raw JSON under `runs/`.

## F1 extremality (R1-R6)

Brute-force oracle over 2^8 = 256 edge subsets, 25 trials per case, n = 7 nodes.
"exact" counts trials where both endpoint gaps are 0 at relative tolerance 1e-12.

| registered prediction | exact trials | max gap_hi | max gap_lo | verdict |
|---|---|---|---|---|
| R1 sum, nonneg W, nonneg X | 25/25 | 0.000e+00 | 0.000e+00 | PASS |
| R2 max, nonneg W, nonneg X | 25/25 | 0.000e+00 | 0.000e+00 | PASS |
| R3 mean | 0/25 | 1.142e-01 | 1.335e-01 | PASS |
| R4 degree-norm | 0/25 | 6.803e-02 | 8.853e-02 | PASS |
| R5 sum, signed W | 11/25 | 1.108e+00 | 8.254e-01 | PASS |
| R6 sum, negative X | 6/25 | 9.195e-01 | 4.907e-01 | PASS |

No falsification criterion triggered. R1 and R2 hold exactly, so the theorem survives.
R3-R6 all break exactness on a majority of trials, so the monotonicity hypotheses have
content rather than being decorative. No trial was voided.

## F2 lattice size (R7)

30 random sequences per length, optional-pair band 0.05 <= p <= 0.90.

| length (nt) | median k | p90 k | max k | median lattice size |
|---|---|---|---|---|
| 40 | 16.0 | 28 | 32 | 10^4.8 |
| 80 | 48.5 | 79 | 95 | 10^14.6 |
| 150 | 86.0 | 135 | 165 | 10^25.9 |

R7 holds. Median k at 150 nt is 86, far above the vacuity threshold of 8. Brute-force
enumeration is out of reach at realistic target lengths, so the two-pass collapse has
practical bite.

## F3 end-to-end certificate

RETIRED. Every magnitude below came from an UNTRAINED, randomly initialised network and
describes the initialisation, not siRNA biology. F3b shows 87.9% of the headline was
relaxation slack, and F6 replaces it with 12.9% using a trained model. Do not quote 110.9%.

20 random 60 nt sequences, untrained monotone model, sum aggregation.
Median relative underestimate of worst-case risk by an MFE-only prediction: 110.9%.
Per-sequence underestimates range from 17.2% to 292.5%. Every certificate cost two
forward passes against lattices of 10^3.3 to 10^17.8 structures.

Caveat carried from `certmp/ensemble.py`: the lattice is a sound outer relaxation of the
Boltzmann ensemble, since not every edge subset is a valid nested secondary structure.
The gap above therefore mixes genuine ensemble uncertainty with relaxation slack. The
conservatism has not yet been quantified. That is the next measurement to make.

## F4 top-k rank stability

Library of 40 candidate sequences at 50 nt, one shared untrained model.
No ranking certified at k = 1, 3, 5 or 10; margins were negative throughout, from
-170.4 at k = 1 to -5024.2 at k = 10.

This is a real negative result, not a code failure. The reachable intervals are wide
enough to straddle every cut, which is what an untrained model over a loose relaxation
should produce. Interval width has to come down, through a trained model or a tighter
ensemble encoding, before rank certification says anything. F4 is not evidence against
the theorem, which F1 tests directly.

## F5 stress audit (NOT pre-registered)

Added after F1 passed, as a post-hoc robustness audit. It can only weaken confidence in
the theorem, never establish it. The registered test remains F1.

An independent standalone reimplementation of the linchpin test, sharing no code with
`certmp`, reproduced F1 digit for digit on all six cases. F5 then randomised everything
F1 held fixed: node count 4-10, k 1-11, 1-4 layers, hidden width 2-23, feature scale
spanning four orders of magnitude, and a random mandatory subgraph.

| aggregation | trials | live | void (dead) | violations | max gap |
|---|---|---|---|---|---|
| sum | 400 | 400 | 0 | 0 | 0.000e+00 |
| max | 400 | 286 | 114 | 0 | 0.000e+00 |

Zero violations in 686 live trials. Endpoint exactness does not depend on any of the
swept quantities.

### The max-aggregation void rate is a real finding

114 of 400 max trials produced an output constant across the entire lattice. Those are
voids, not confirmations, so they were excluded. Chasing the cause: the void rate is
independent of k, and driven entirely by the density of the mandatory subgraph. At
n = 8, k = 5:

| mandatory edges | sum live | max live |
|---|---|---|
| 0 | 1.00 | 1.00 |
| 4 | 1.00 | 1.00 |
| 8 | 1.00 | 0.97 |
| 12 | 1.00 | 0.73 |
| 16 | 1.00 | 0.33 |
| 20 | 1.00 | 0.05 |
| 23 | 1.00 | 0.00 |

A dense mandatory core already saturates the row-wise maximum, so optional edges cannot
move it. Max aggregation is therefore exactly certifiable but frequently uninformative:
the certificate collapses to worst case equals best case. Sum never degrades this way.
This is a limitation of S2 that the pre-registration did not anticipate, and it belongs
in any write-up of the aggregation asymmetry.

### The failure mode does not reach the RNA application

On the 20 folded 60 nt sequences from F3, mandatory pairs are those with p > 0.90. A
base pairs with at most one partner, so the mandatory subgraph is a matching: maximum
node degree 1, median density 0.28%, maximum 0.68%. That is far below the roughly 40%
density where max aggregation begins to saturate. Max stays live on real ensembles.

## F3b relaxation slack (follow-up to F3, added 2026-09-04)

F3's headline number was not what it appeared to be. The edge lattice is an outer
relaxation, so its maximum can be attained at an edge subset that no real secondary
structure realises. This experiment separates genuine ensemble uncertainty from
relaxation slack by drawing 1000 Boltzmann samples per sequence with `RNA.pbacktrack`,
seeded for reproducibility, on the same 20 sequences and the same models F3 used.

### Most of F3's reported gap is slack, not uncertainty

| quantity | median |
|---|---|
| F3 reported gap, MFE vs lattice worst case | 110.9% |
| honest gap, MFE vs sampled ensemble maximum | 13.4% |
| share of F3's reported gap that is pure relaxation slack | 87.9% |

The lattice maximum runs a median 1.82x the true sampled ensemble maximum, ranging from
1.10x to 3.23x. The real effect of ensemble uncertainty on this model is 13.4%, not
110.9%. F3's number should not be quoted.

### The upper bound is sound in scope but narrower than claimed

Only a median 73.8% of sampled structures lie inside the lattice, worst sequence 43.1%.
The band keeps pairs with probability at least 0.05 and discards the rest, so a sampled
structure using a discarded pair is not covered by the certificate at all. The bound held
on 100% of samples in this run regardless, but that is an empirical observation and not a
guarantee. Building the lattice with the lower band at zero would make it a guarantee, at
the cost of a larger k.

### The lower bound is not sound, and that is a defect

The lattice minimum forces every mandatory pair present. Real structures do not: a median
96.6% of samples contained all mandatory pairs, worst 83.0%. The certified minimum was
violated by an actual sampled structure on 2 of 20 sequences. Lowering the band does not
fix this, because the problem is the forcing, not the discarding. The lattice minimum is
a lattice diagnostic, not an ensemble bound.

Docstrings in `certmp/ensemble.py` and `certmp/certify.py` claimed soundness in both
directions. Both have been corrected to match what was measured, as has README.md.

### The MFE structure is always inside the lattice

Every minimum-free-energy base pair fell inside mandatory union optional on all 20 of 20
sequences, with zero failures. That check passes cleanly.

## F6 trained model (replaces every magnitude in F3 and F3b)

F3 and F3b used a randomly initialised network. This trains one on real measured siRNA
efficacy and repeats the comparison.

### Dataset

The Huesken et al. 2005 table was obtained from the DSIR redistribution (Vert et al.,
BMC Bioinformatics 7:520, 2006) and cross-checked against two independent redistributions,
which agree exactly on both sequences and values. It passes every integrity check in
`certmp.data.verify` unaided: exactly 2431 rows, exactly the published 2182/249 split,
alphabet ACGU, zero duplicate sequences, zero train/test overlap.

Two pre-registered expectations were wrong and are recorded as amendments in
`certmp/data.py` rather than quietly edited. Sequence length is 21, not 19, and the final
two nucleotides vary across rows so they are measured sequence rather than a constant
overhang. Efficacy reaches 1.341, above the pre-registered ceiling of 1.2, because the
published values are normalised inhibition and legitimately exceed 1.0. The integrity
check for the efficacy range had never been implemented; it is now.

### Held-out performance

| split | n | Spearman | Pearson |
|---|---|---|---|
| published 2182/249 | 249 | 0.592 | 0.556 |
| gene-disjoint | 290 | 0.543 | 0.531 |
| Huesken 2005 published baseline | 249 | not reported | 0.66 |

The published split shares all 30 annotated genes between train and test, so it measures
generalisation to new siRNAs against known targets. A gene-disjoint split, the harder and
more relevant question for a design tool, costs only about 0.05 Spearman, so the leakage
is not inflating the headline much.

The shortfall against the published 0.66 is the price of certifiability. The baseline was
a feedforward ensemble on hand-built features with no sign constraints. A monotone
sign-constrained GNN is a strictly smaller hypothesis class, and it recovers most but not
all of the signal.

Two defects had to be fixed before this number meant anything, and the first attempt
produced a Spearman of -0.095. The readout was six orders of magnitude off the target
scale, so training spent every epoch correcting scale rather than learning; the model now
carries a monotone affine readout initialised from training statistics. Separately, sum
aggregation is permutation invariant, so plain one-hot nucleotide features made the model
a bag of nucleotides that could not express position-specific preferences, which is most
of the siRNA signal. Features are now one-hot over position and nucleotide jointly, which
stays non-negative and satisfies H4. Both the non-negative scale and the shift are
order-preserving, and `tests/test_sanity.py` now checks that they leave extremality and
torch/numpy parity intact.

### The trained comparison

| quantity | untrained (F3/F3b) | trained (F6) |
|---|---|---|
| MFE vs lattice worst case | 110.9% | 12.9% |
| MFE vs sampled ensemble maximum | 13.4% | 11.4% |
| relaxation slack ratio | 1.82x | 1.00x |

F3's headline 110.9% was an artifact of random initialisation and is now retired. With a
trained model the MFE-only prediction understates worst-case risk by about 13%, and
almost all of that is genuine ensemble uncertainty rather than relaxation slack.

The two slack ratios are NOT comparable and must not be quoted side by side as a
reduction. F3b folded 60 nt sequences with a median of 30 optional pairs; F6 folds
21 nt siRNAs with a median of 6. The lattice is nearly tight at k = 6 simply because
there is almost nothing to relax.

### Limitation: the certification story does not live on the siRNA itself

An isolated 21 nt siRNA has a median of 6 optional pairs, so the uncertainty set holds
about 64 structures and brute force would be entirely affordable. Three of the held-out
sequences scanned had no optional pairs at all. The exponential-to-constant collapse that
F2 established needs the longer mRNA target context, where k reaches 86 at 150 nt, not
the siRNA duplex. F6 demonstrates that a certifiable model can be trained to useful
accuracy; it does not demonstrate that certification is needed at this sequence length.

## F7 sound certificate and the price of soundness (step 1)

f3b found the certificate was unsound in both directions. This fixes it and measures what
the fix costs. Model is the length-generalizing monotone network trained on real Huesken
efficacy, held-out Spearman 0.385, not an untrained network. 20 sequences at 60 nt, 1000
seeded Boltzmann samples each, 20000 sampled structures in total.

The fix is structural. Setting mandatory to empty means nothing is forced present, so the
empty endpoint is a subset of every real structure. Setting the probability floor to zero
discards no pair, so the lattice is the full power set of base pairs and every secondary
structure is inside it. Monotonicity then makes the upper endpoint a genuine bound on the
whole ensemble.

### Coverage versus tightness

| floor | canonical only | median k | coverage | worst coverage | lattice max / ensemble max | bound violations |
|---|---|---|---|---|---|---|
| 0 | no | 1770 | 100.0% | 100.0% | 179.08 | 0 / 20000 |
| 0 | yes | 610 | 100.0% | 100.0% | 29.25 | 0 / 20000 |
| 1e-6 | yes | 414 | 100.0% | 100.0% | 15.14 | 0 / 20000 |
| 1e-4 | yes | 190 | 99.7% | 99.3% | 4.69 | 0 / 20000 |
| 1e-3 | yes | 120 | 98.2% | 96.6% | 2.57 | 0 / 20000 |
| 1e-2 | yes | 61 | 91.0% | 79.2% | 1.48 | 0 / 20000 |
| 0.05 | yes | 34 | 73.8% | 43.1% | 1.13 | 15 / 20000 |
| 0.10 | yes | 26 | 59.5% | 36.3% | 1.01 | 199 / 20000 |
| 0.20 | yes | 20 | 41.9% | 8.1% | 0.98 | 1017 / 20000 |

Restricting the lattice to the six canonical pairs with the minimum hairpin loop costs
nothing in soundness, because ViennaRNA's model cannot form any other pair, and it cuts
the lattice to a third of the edges and the bound from 179x to 29x. Only the backbone is
mandatory; it is not a base pair, every structure has it, and the model saw it in training.

The curve is written to `coverage_tightness.png` in the run directory.

At floor zero the bound holds on 100% of sampled structures with zero violations, as the
theorem requires. The old default of 0.05 covered only 73.8% and was not merely
theoretically unsound: 15 real sampled structures scored above its certified maximum. At
0.20 more than a thousand did.

### Soundness costs a factor of 29, and that is the headline

The sound bound runs a median 29x the true sampled ensemble maximum. A certificate that
passes at floor zero is therefore strong evidence, but one that fails may mean only that
the relaxation is loose. The tightness the earlier experiments reported, 1.13x at floor
0.05, was bought by excluding a quarter of the ensemble from the claim.

The middle of the curve is where a practitioner would actually sit. A floor of 1e-3 covers
98.2% of sampled structures at 2.57x, and 1e-2 covers 91.0% at 1.48x. Neither is sound.
There is no floor in this data that is both sound and tight, and that is the real finding:
the edge lattice is too coarse a relaxation of a nested secondary structure to be both.

### The lower bound is gone

`certify_threshold` no longer returns a best case, and `certify_topk` is withdrawn and
raises. Exact top-k rank certification needs a sound lower bound per candidate, and certmp
no longer claims one, so F4 is retired rather than left standing on an invalid bound. It
certified nothing when it ran. F3 is retired too: untrained model, unsound banded lattice.
`reach.exact_interval` still returns both lattice endpoints, because the extremality
theorem is a statement about both and F1 and F5 test both. That is lattice arithmetic, not
an ensemble claim.

## F8 matched baseline, the cost of certifiability (step 2)

The baseline differs from the monotone model in exactly one respect, the sign constraint.
Identical features, published split, validation set, optimiser, learning rate, batch size,
epoch budget, early stopping and seeds, all through the same `certmp.train.fit`. Five seeds
per configuration, reported as mean and standard deviation.

| features | model | Spearman | Pearson |
|---|---|---|---|
| positional | monotone (certifiable) | 0.592 ± 0.004 | 0.557 ± 0.004 |
| positional | unconstrained | 0.594 ± 0.013 | 0.588 ± 0.013 |
| generic | monotone (certifiable) | 0.389 ± 0.012 | 0.372 ± 0.010 |
| generic | unconstrained | 0.603 ± 0.007 | 0.594 ± 0.008 |
| positional, gene-disjoint | monotone (certifiable) | 0.544 ± 0.002 | 0.531 ± 0.001 |
| positional, gene-disjoint | unconstrained | 0.583 ± 0.023 | 0.570 ± 0.024 |

### Certifiability is nearly free, but only with an expressive non-negative encoding

With one-hot (position, nucleotide) features the sign constraint costs 0.002 Spearman,
which is inside the seed-to-seed spread, and 0.031 Pearson. On a gene-disjoint split it
costs 0.039 Spearman. That is the headline: the extremality theorem is close to free here.

With the compressed generic encoding it costs 0.214 Spearman, more than a third of the
achievable signal. The unconstrained model reaches 0.603 on those same features, so the
signal is present and it is the non-negativity that cannot extract it.

The mechanism is straightforward. Under one-hot (position, nucleotide), every combination
gets its own non-negative weight, so a preference for one nucleotide at a position is
expressed by giving the alternatives smaller weights, and no negative weight is needed.
Under a compressed encoding, features share weights and the model needs subtraction it is
not allowed to perform.

The practical consequence is a design rule rather than a limitation: pay for
certifiability in the feature encoding, not in accuracy. It also qualifies F7 and F9,
which use the generic encoding out of necessity because they span several sequence
lengths, and therefore run on a model handicapped by roughly 0.21 Spearman.

## F9 sound certificate on real mRNA target sites (step 3)

f6 showed the certificate was being run in the wrong place. An isolated 21 nt siRNA has a
median of 6 optional pairs, so brute force is affordable and the theorem earns nothing.
This runs it on the mRNA target site instead, extracted from GENCODE release 47 transcripts
for the Huesken genes.

GENCODE matched 28 of the 30 annotated gene symbols across 219 transcripts, after mapping
four legacy symbols to current names. 1720 of the 2431 siRNA target sites were located
exactly inside a transcript. 223 failed because the gene is rat and absent from a human
build, and 488 because the site did not appear in any GENCODE transcript of that gene. 24
sites were sampled for the certificate, with 300 Boltzmann samples per window.

| window | median k | lattice size | MFE point estimate | sampled ensemble max | sound worst case | slack | coverage |
|---|---|---|---|---|---|---|---|
| 50 nt | 408 | 10^123 | 1.0 | 1.1 | 23.7 | 21.5x | 100% |
| 100 nt | 1744 | 10^525 | 1.5 | 1.5 | 146.3 | 94.3x | 100% |
| 150 nt | 3932 | 10^1183 | 1.9 | 2.0 | 446.4 | 222.8x | 100% |

Every certificate cost two forward passes. The bound was sound on every site at every
length, with 100% of sampled structures inside the lattice.

### The collapse is real here, and it is enormous

At 150 nt the uncertainty set holds about 10^1183 edge subsets and two forward passes
bound all of them. This is the regime the theorem was for. F2 predicted it from base-pair
counts; F9 confirms it on real target sites with a trained model and a sound lattice.

### But soundness and usefulness move in opposite directions with length

The slack grows faster than the window: 21.5x at 50 nt, 94.3x at 100 nt, 222.8x at 150 nt.
The certified worst case at 150 nt is over 200 times the largest value any of 300 sampled
structures actually produced. A threshold certificate that passes at that looseness is
still meaningful, since the bound is genuine, but almost nothing will pass it.

Meanwhile the quantity the certificate is meant to protect against stays small. The MFE
point estimate understates the sampled ensemble maximum by only 4.4% to 5.1% across all
three lengths. On this model and these targets, ensemble uncertainty is a few percent while
the sound bound overshoots by two orders of magnitude.

That is the central tension the theorem now faces, stated plainly. The exponential-to-
constant collapse is real, dramatic, and lands exactly where the biology needs it. The
relaxation that makes it sound is too coarse for the bound to be actionable at that scale.
Closing that gap needs a tighter sound relaxation of nested secondary structure, not a
better network. Nesting and the one-partner-per-base constraint are exactly the structure
the edge lattice throws away.

## Status

Theorem survives its falsification attempt, and survives a randomised robustness audit
it was not pre-registered against.

The model is trained and real: Spearman 0.592 on the published split, 0.543 gene-disjoint,
against a published Pearson baseline of 0.66. F3's 110.9% is retired as an initialisation
artifact and replaced by 12.9%.

The certificate is now sound, one-sided, and running in the right place. F9 bounds 10^1183
structures in two forward passes on real GENCODE target sites, with every sampled structure
inside the lattice.

What remains open is tightness, not soundness. The sound bound overshoots the sampled
ensemble maximum by 29x at 60 nt and 223x at 150 nt, while the effect it guards against is
about 5%. The edge lattice discards nesting and the one-partner constraint, and that is
where the looseness comes from. A tighter sound relaxation is the next piece of work, and
it is a modelling problem rather than a network problem.
