# Pre-registration -- certmp
Frozen before results. Append AMENDMENTS with date/justification; never edit above.

## Claim
For a message-passing GNN with (i) non-negative message/update weights and biases,
(ii) monotone non-decreasing activations, (iii) non-negative node features, and
(iv) non-negative readout, the exact reachable output interval over an uncertainty set
that is a LATTICE of edge subsets (mandatory edges always present, optional edges free)
is attained at the lattice endpoints, and is therefore computable in TWO forward passes
rather than 2^k.

## The aggregation asymmetry (the actual contribution)
| aggregation | feature-box perturbation | edge-lattice perturbation |
|-------------|--------------------------|---------------------------|
| sum         | exact                    | exact                     |
| max         | exact                    | exact                     |
| mean        | exact                    | NOT exact                 |
| degree-norm | exact                    | NOT exact                 |

Mean and symmetric degree normalisation divide by a function of the DEGREE, which the
edge perturbation itself changes. Hence GCN-style normalisation forfeits structural
certification while retaining feature certification.

## Registered predictions
- R1 sum + nonneg W + nonneg X: relative endpoint gap == 0 on every trial.
- R2 max + nonneg W + nonneg X: relative endpoint gap == 0 on every trial.
- R3 mean: gap > 0 on a majority of trials.
- R4 degree-norm: gap > 0 on a majority of trials.
- R5 signed weights: gap > 0 (monotonicity is what buys exactness).
- R6 negative node features: gap > 0 (the X >= 0 hypothesis is necessary, not cosmetic).
- R7 On real ViennaRNA Boltzmann ensembles, the optional-pair count k grows with target
  length such that 2^k is intractable at realistic lengths (>= 50 nt), establishing that
  the exponential-to-constant collapse has practical bite rather than being vacuous.

## Falsification criteria
- If R1 or R2 fails (a strict interior optimum exists), the theorem is FALSE. STOP.
- If R7 fails (k stays small, say median k < 8 at 150 nt), the result is TRUE BUT VACUOUS,
  because brute force would already be cheap. Reframe or abandon.
- If R3-R6 all PASS extremality (i.e. everything is exact), then monotonicity is not the
  operative condition and the characterisation has no content. STOP.

## Void conditions
A trial is VOID, not a confirmation, if: the model output is constant across the whole
lattice (dead network); the optional-edge set is empty (k = 0, nothing to certify);
any non-finite value appears; or gaps are compared on an absolute rather than relative scale.

## Status
UNRUN.

## AMENDMENTS
2026-09-04: first full run executed. Nothing above this line was edited.
Results recorded in RESULTS.md and runs/. Makefile gained `export PYTHONPATH := .`
because `python3 tests/test_sanity.py` put tests/ rather than the repo root on sys.path;
this is an invocation fix and changes no experimental logic.
