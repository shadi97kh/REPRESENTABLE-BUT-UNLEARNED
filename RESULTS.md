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

## Status

Theorem survives its falsification attempt. The application is not yet certified.
