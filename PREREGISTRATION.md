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
2026-09-04: added experiments/f5_stress_extremality.py, a POST-HOC robustness audit of
the two surviving cases S1/S2. Not a registered prediction and not evidence for the
theorem; it exists to falsify, and it did not. It also surfaced an unanticipated
limitation of max aggregation (void under a dense mandatory subgraph) recorded in
RESULTS.md. No experiment above this line was modified.
2026-09-04: two EXPECTED values for the Huesken dataset in certmp/data.py were wrong on
first contact with the real distribution and are amended there with justification, not
silently edited. Sequence length 19 -> 21, because the distributed sequences are 21 nt and
the final two nucleotides vary. Efficacy ceiling 1.2 -> 1.4, because the observed maximum
is 1.341 and the published values are normalised inhibition that exceeds 1.0. Every
diagnostic check passed unaided: 2431 rows, the 2182/249 published split, ACGU, no
duplicates, no train/test overlap. The efficacy-range check had never been implemented and
now is. Added experiments/f3b_relaxation_slack.py and experiments/f6_trained.py, neither
pre-registered. F3's headline number is retired as an initialisation artifact; see
RESULTS.md.

2026-09-04 (corrections round). Nothing above this line was edited. All items POST-HOC.

C1 BOUNDED DEGREE REMOVED. An earlier draft planned to assume bounded degree to avoid
   colliding with Saelzer & Lange (ICLR 2023). That assumption is unnecessary: their
   undecidability quantifies over unbounded graph FAMILIES, whereas this theorem is over
   a fixed finite graph with a finite optional-edge set, where reachability is decidable
   by enumeration. Bounded degree belongs in related work, not in the hypotheses.

C2 COMPLEXITY CLAIM NARROWED. The defensible statement is exact reachability in two
   network evaluations, independent of k, contrasted with 2^k enumeration. Claims
   resting on the NP-hardness of third-party verifiers are not load-bearing.

C3 AGGREGATOR CHARACTERISATION GENERALISED (c3_aggregator_characterization.py).
   Endpoint exactness holds iff the aggregator is monotone under multiset inclusion, in
   EITHER direction. Antitone aggregators (min) are certifiable with the endpoint roles
   swapped. reach.py now dispatches on model.endpoint_direction(). This strictly
   generalises the original four-row sum/max/mean/degnorm table and admits logsumexp
   and min. Verified on 7 aggregators, including the direction of the maximum.

C4 AFFINITY MEASURED (c4_affinity.py). Under H1-H5 with ReLU every preactivation is
   non-negative, so ReLU never clips and a linear-aggregator network is exactly affine
   (additivity violation 2e-16). This is an expressivity limitation of the ReLU
   instantiation, not of the certified class; tanh and sigmoid satisfy H3 and restore
   nonlinearity. Any near-free-certifiability claim must carry the qualifier that the
   task is near-linear in the encoding used. Every model in f6 through f10 was affine.

C5 SLACK SCALING EXPLAINED (c5_slack_scaling.py). The exponent is set by network DEPTH,
   not by RNA. On a synthetic matching family containing no RNA the fitted exponent is
   0.912 per layer at depths 1 to 4 with R^2 0.9996. The 1.93 exponent measured on RNA
   with a 2-layer network matches 2 x 0.912 = 1.824 and was therefore an architectural
   property misread as thermodynamics. Slack is also unbounded in problem size, and is
   exactly 1.0 on join-closed families, verified by brute force rather than asserted.

C6 CLOSURE (c6_closure.py). Endpoint exactness survives non-negative residual
   connections, depths 1 to 5, and every monotone activation: 52/52 live configurations,
   0 violations, 2 voids reported separately. A non-monotone activation breaks it in all
   4 control configurations, so the test discriminates.

Experiment files are named by correction ID (c3..c6) because f7 through f10 were already
taken by the soundness, baseline, target-context and maximal-element experiments.

2026-09-04 (kappa round, additive). POST-HOC. Adds certmp/kappa.py and c7 to c9. No
existing file was overwritten.

C7 The relaxation gap of the certified class equals kappa_L = W_L(top)/max_F W_L, a ratio
   of length-L walk counts, exactly in the no-bias uniform-feature regime (max relative
   error 2.12e-16) and as a tight upper bound with biases. Join-closed families give
   exactly 1.0. Verified by exhaustive enumeration, not sampling.
C8 On the same GENCODE windows as f9/f10, kappa_2 bounds the measured slack at every
   length and reproduces its growth (4.17x, 2.22x against 4.39x, 2.36x measured), while
   overshooting the level by 5.2x through feature and bias dilution. kappa_L = kappa_1^L
   to 4.6%, so C5's depth proportionality is combinatorial, not empirical.
C9 Equality holds only in the affine uniform-feature regime; kappa_L still upper-bounds
   the gap under varied features, max aggregation and tanh.

Two fixes were required to make c8 run and be reproducible, both recorded here rather than
applied silently. walk_count now uses L matrix-vector products instead of a matrix power,
mathematically identical and verified against the reference, without which c8 would not
finish. RNA.cvar.rand_seed does not exist in this ViennaRNA binding and was silently
swallowed by a try/except, leaving sampling unseeded; it is replaced by RNA.init_rand.
data/gencode_windows.txt is generated from the same sites f9 and f10 used, so c8 ran on
real windows rather than its random fallback.
