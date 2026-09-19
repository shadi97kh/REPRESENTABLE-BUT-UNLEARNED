> Package note — copied historical evidence. Mathematical prose and formulas are retained; links to excluded repository artifacts are rendered as inactive references, and any link changes are recorded in the package manifest. This is not a new proof or an edit to the original file.

**Decision: incremental-known. A larger campaign with the present method is not scientifically justified.**

This review used the completed artifacts only. It inspected the learner, evaluator, correctness-test source and saved test output, mathematical notes, primary-source comparisons, frozen configuration, saved predictions and fold records, raw results, all three existing PNG figures, and resource/preservation records. Two independent static agent reviews covered mathematics and experiment design. No fits, numerical experiments, test execution, model evaluation, resampling, or figure regeneration were performed. Hashing and arithmetic summaries of existing JSON/CSV records were bookkeeping. The original code, reports, predictions, figures and ledgers were preserved; this is a new audit note.

**Strongest defensible contribution and closest equivalence.**

The work provides a carefully specified five-source-law interaction problem, globally admissible fitted profile families, a correct transformation-law score and target derivative, an implementable source-only correction, and a reproducible pilot with useful negative evidence. The score at fixed outcome is

`s_v(F(z)) = z*u_v(z) - u_v'(z)`, with `u_v = dot F_v/F'`.

This model-specific calculation, its Gaussian integration-by-parts justification, the ideal conditional MSE identity, and the admissible finite-path derivative survive scrutiny. They specialize established likelihood/inverse-problem machinery. The finite correction solves a Tikhonov-regularized score-Riesz problem; the next direction maximizes a stabilized generalized Rayleigh quotient and is appended by Galerkin enrichment. No distinct learned finite-sample rate has been proved.

This is an algebraic equivalence after freezing the score dictionary and metric, not a claim that native regression RieszNet can be run unchanged on these observations. [RieszNet/ForestRiesz](https://arxiv.org/html/2110.03031v3) already learns representers for orthogonal functional correction. [Trabs's generalized-score theory](https://arxiv.org/abs/1307.6610v2) supplies the adjoint-range formulation for inverse statistical experiments. [Strongly identified functionals of weakly identified functions](https://arxiv.org/html/2208.08291v3) already separates functional recovery from nuisance recovery. The potentially distinct task is to establish and attain an original-model nonlinear contrast rate under the particular five-law coverage and weak separation, with learned geometry and honest approximation control. The current work does not complete that task.

[DExtrI](https://arxiv.org/html/2608.19849v1) is also close in shared-latent profile fitting with proper scores. Its axis-interval observation model and identification assumptions differ from the stipulated known warp, two private unit anchors and five discrete laws here. Those differences require a new proof; they do not by themselves constitute one. Its broader appendix rigidity results must be considered, rather than relying only on a main growth assumption. The source matrix's qualifications about native IV, bilinear-completion and PDE comparators are appropriate: their observations or operator assumptions are unavailable here, while their applicable spectral and residual algebra remains prior art. This was a bounded source audit, not proof that every potentially relevant paper has been excluded.

**1. The learner does not currently satisfy everything needed for its claimed guarantees.**

The projected parameter constraints and analytic profile derivative envelopes are valid in real arithmetic. The correction weights and independent evaluation role support the stated conditional-risk calculation for a well-defined estimator. These are meaningful positive checks.

There is, however, a concrete mismatch between the executable and the coarse global finite-risk statement. profiles.py (omitted repository reference: `../../../interaction_recoverability/learned_target/profiles.py#L111`) fixes the inverse bracket to `[-24,24]` and raises `ArithmeticError` outside the fitted response image of that interval. correction.py (omitted repository reference: `../../../interaction_recoverability/learned_target/correction.py#L23`) has no fallback. The admissible Gaussian-latent source law has support `(0,infinity)`. Conditional on any fitted model, outcomes outside that finite image have positive probability. Therefore the literal executable is not a total estimator on the statistical sample space, whereas theorem.md (omitted repository reference: `theorem.md#L53`) uses a globally defined ideal inverse to bound its MSE and calls the result an "implemented-form" bound.

The envelope argument can support the ideal total estimator. Applying it to the executable requires an expanding inverse bracket or an explicit fallback, together with an error/failure contribution. Conditioning on successful numerical evaluation alone does not supply the unconditional guarantee. No saved dataset encountered this failure, so this objection does not invalidate its recorded predictions. The residual check `abs(F(z)-y)/(1+y)` also does not certify the log-response-to-index tolerance invoked in the proof, particularly near zero. The simplified variance term `C(V)/N_E` further presumes the sampling weights used by the implementation; arbitrary weights require the preceding sum of `pi_a^2/n_a` terms.

The sharper guarantee is explicitly open for deeper reasons: no shrinking source-supported radius, learned score-operator transfer bound, original-class direction-approximation bound, uniform nonlinear-bias bound, or tuning oracle inequality is established. There is a specific approximation obstruction to an unstated neural-consistency assumption: every allowed tanh network has a convergent `g(x)/exp(x)` as `x` tends to positive infinity, whereas the oscillatory truth retains sine modulation. The truth in evaluate.py (omitted repository reference: `../../../interaction_recoverability/learned_target/evaluate.py#L17`) is therefore outside this fixed neural family. A target may still be recoverable despite nuisance misspecification, but that requires a target-specific proof. Fixed-node CRPS training likewise does not supply a demonstrated rate from the fitted discrete objective to the continuous true-law geometry.

**2. The direction search extends beyond network parameter derivatives, but only within a fixed finite bank.**

geometry.py (omitted repository reference: `../../../interaction_recoverability/learned_target/geometry.py`) builds 36 physical directions: four coefficient directions and 16 analytic profile perturbations for each component. The sine, Gaussian and tanh perturbations are not restricted to differentiating the fitted four-unit network's parameters. The coupled finite paths can therefore leave its architecture while remaining inside the larger admissible profile class for the checked step sizes.

All adaptive directions nevertheless lie in this same 36-dimensional span. The fixed `neural_rich` comparator already spans it. Search does not examine arbitrary analytic profile perturbations, globally different source-compatible profiles, alternative latent laws or unknown warps. An exact optimum for the numerical finite-bank quotient is not a whole-class uncertainty characterization.

There is an additional interpretation problem. For a direction `v` in the retained correction span `V`, the regularized normal equations imply

`v^T(b-Gc) = lambda * <v,c>_M`.

Thus the residual generally remains nonzero inside the current span because of deliberate ridge bias. enrichment.py (omitted repository reference: `../../../interaction_recoverability/learned_target/enrichment.py#L13`) searches the full residual without removing this component. Its witness is a valid regularized residual diagnostic, but it does not isolate omitted uncertainty. A large value need not mean the learner overlooked a new direction. This also helps explain why interpretation requires more than a nonzero witness.

**3. Finite witnesses were not presented as certified upper bounds.**

The proof, code labels and figures correctly distinguish finite witnesses from full-class bounds. Of 144 evaluated paths, 124 passed the reported source-compatibility check; this is not coverage of a confidence set. The basis endpoint search returns approximate inner witnesses in a restricted family, not an original-class confidence interval or certified extrema. Failure to find a larger contrast does not bound the remaining contrast range.

The finite-path diagnostic uses different quadrature resolution from the operational target/geometry calculation. Its reported nonlinear remainder consequently includes numerical baseline differences as well as path curvature. The notes acknowledge this. It should not be read as an isolated, uniform bound on `B''`.

**4. Oracle information is excluded from fitting, but several quantities needed by a useful theorem remain unknown.**

The learner receives only the five source outcome arrays and stipulated public model information, including the Gaussian latent law, known warp/separation and anchor design. No generator latent observations, true profiles, true shifts, held-out joint outcomes or true target values enter fitting or selection. Public integration nodes are not observed latent variables. Configuration/code commitments precede final generation, and predictions are hashed before evaluation-only truth access. This is an auditable interface separation, not operating-system isolation from a malicious process.

The `.02` term in selection_score (omitted repository reference: `../../../interaction_recoverability/learned_target/correction.py#L37`) is a heuristic scale, not an estimated nuisance radius. The finite-linear spectral calculation assumes fixed geometry, fixed lambda and a supplied radius. Turning it into a useful learned guarantee additionally needs true-to-fitted operator/variance transfer, omitted-direction control, selection error control, and a uniform nonlinear remainder. The huge global envelope bound does not make these terms shrink. Small observed remainders do not replace those assumptions.

**5. Correction algebra and sample separation are substantially correct.**

The fixed-outcome score includes the inverse-index derivative and the `F''*dot F/(F')^2` term. Target differentiation includes both profiles and all four shifts. The plus sign in `plugin + correction` is correct. Evaluation weights are `pi_a=n_a,E/N_E`; the variance calculation uses `sum_a pi_a^2 * sample_variance_a/n_a,E`. Scores are centered separately under the fitted source laws, using the declared numerical integration rule. Centering error is retained as bias rather than silently treated as exact.

Each rotation has disjoint fit, selection and evaluation observations. Checkpoints, directions and lambda are chosen without that rotation's evaluation outcomes. The three cyclic estimates are dependent, and the proof correctly uses a Jensen bound instead of claiming three independent replicates or a calibrated confidence interval.

One generic extension caveat: selection-score averaging defaults to the selection sample proportions while the Gram matrix uses evaluation proportions. They coincide for the actual balanced per-stratum allocation, so this is not a demonstrated bug in these 24 datasets. Unequal-stratum extensions need consistent weighting or a revised derivation.

**6. Information was comparable; capacity, optimization and cost comparisons have qualifications.**

The strongest same-fit controls share nuisance fits, source observations, geometry and lambda grid. In fact, adaptive and rich selected the same lambda in all 72 rotations. This is a useful control against attributing their difference to a different selected ridge value.

The random control is not fully matched in effective rank. Candidate acceptance uses a residual M-norm threshold of `1e-7`, but the later Gram eigendecomposition discards directions below a relative eigenvalue threshold of `1e-9`. An accepted tiny addition can disappear. All 288 adaptive appends passed the precheck, yet adaptive ended at rank 12 in 37 rotations, 11 in 32, and 10 in 3. Random ended at rank 12 in all 72. Thus the adaptive/random comparison also changes effective dimension and implicit truncation; equal append rounds do not isolate targeted direction choice. The earlier report disclosed adaptive rank loss, but did not make this random-control asymmetry explicit.

The pooled controls have one fit and use the final checkpoint, while shared/corrected methods average three separately initialized, source-checkpoint-selected fits. Their gap confounds correction with ensembling, sample roles and checkpoint selection. The same-fit plug-in and full-bank controls are the more informative ablations.

The direct-target baseline received the same source information and more recorded CPU, but its finite-basis, penalty-based endpoint search is not a strong solution of the full confidence-set problem. The best feasible point was at the step-40 boundary for 30 of 48 endpoints. This does not prove failed convergence, but it rules out treating feasibility or job completion as an optimization certificate.

The recorded per-method CPU bars are component costs, not complete end-to-end runtimes. experiment.py (omitted repository reference: `../../../interaction_recoverability/learned_target/experiment.py#L77`) adds selection/construction time but does not assign the subsequent adaptive finite-path diagnostics to its method total. Target integration and applying the final correction also fall outside these method counters. The dataset totals and resource ledger include these operations; this finding does not establish an allowance-accounting shortfall. The existing description "fully charged" should be read only as fully charging the measured shared components, not every operation. Every corrected comparator computes the complete 36-dimensional geometry first, including the fixed eight-direction method, so no avoided-full-bank computational advantage is demonstrated.

Unfavorable outcomes were retained. The saved final experiment contains all 24 datasets, all ten methods and all 240 prediction rows, with no failed final fit hidden from the table. Development/warning records and earlier negative work remain. The only nonzero wrapped job was the intentional watchdog test. No post-evaluation retuning or rerun was found in the ledger.

**7. The observed improvement does not isolate the proposed mechanism.**

The existing raw rows give the following descriptive results; these are not new experiments or significance tests.

| Existing-artifact comparison | Result |
|---|---:|
| Adaptive aggregate MSE | 0.4793938283 |
| Full-bank neural aggregate MSE | 0.4970740215 |
| Relative reduction versus full bank | 3.5569% |
| Frozen descriptive screening threshold | 10% |
| Adaptive lower squared error versus full bank | 13/24 datasets |
| Adaptive lower squared error versus fixed / random | 9/24 / 10/24 datasets |
| Adaptive/rich identical selected lambda | 72/72 rotations |
| Median same-direction witness after/before | 0.9985769849 |

Two datasets, `d012` and `d013` (localized truth, delta .1, n=256), account for approximately 96% of the net paired squared-error gain over the full-bank control: their combined reduction is about .407805 out of .424325 across all 24 datasets. This is a concentration statement about the saved results, not a post hoc reason to remove those datasets. Three replicates per cell, a small aggregate gain, rank mismatch, heuristic selection and a residual that mixes omitted directions with ridge bias do not establish that uncertainty-aware direction discovery caused a robust advantage. Passing the 10% threshold alone would not have proved novelty either.

All three original figures were inspected. The first retains individual replicates; the second supports same-fit comparisons but needs the CPU qualification above; the third visibly shows many witness pairs near the diagonal. No figure supplies an uncertainty certificate.

Other retained diagnostic values are:

| Diagnostic | Saved value / interpretation |
|---|---|
| M refinement spectral-norm difference, 4096 vs 2048 nodes | 1.368362873580413; not a certified Sobolev-norm approximation |
| Maximum score-Gram refinement difference | 0.0010506066309624639 |
| Maximum fitted-score centering magnitude | 0.000015847774228128836 |
| Maximum selected linear-system relative residual | 1.2345301667400524e-15; certifies the solved numerical system only |
| Maximum selected regularized condition number | 101.09694066775846 |
| Maximum absolute finite-path B | 0.10207009503408186 |
| Maximum reported nonlinear remainder | 0.0012669990194164732; includes numerical baseline effects |

**8. Principal unresolved objection and scientific next gate.**

The central gap is an attainable finite-sample bound for the target over the complete five-source-law experiment with unknown profiles. The finite fitted residual cannot supply it without an honest bound on uncertainty outside the bank and on nonlinear source-to-target transfer. This objection remains even after repairing the inverse and rank thresholds. The retained known-profile complete-experiment lower bound is useful subclass evidence; no matching constructive unknown-profile upper rate was supplied.

A substantive contribution would need either a model-specific target modulus/rate and an estimator attaining it under verified conditions, or a clearly distinct operation with persuasive evidence against equally informed established alternatives. Before another empirical campaign, the estimator must be defined on the stated sample space; numerical and approximation assumptions must be explicit; ridge residual and omitted-direction effects must be separated; effective-rank and checkpoint/ensemble controls must be matched; and complete cost attribution and meaningful endpoint optimization checks must be specified. These are prerequisites, not authorization or a recommendation to run training now. More repetitions of the current implementation would not resolve the main objection.

**9. Nothing here establishes a biological siRNA interaction.**

The siRNA mapping (omitted repository reference: `sirna_mapping.md`) correctly requires a measured reference/single/single/joint quartet on the same sequence and assay background, with exact guide/passenger identities, chemistry and positions, linkages/stereochemistry, dose, time, cells, batches, shared controls and biological replicate provenance. Joint outcomes would need independent held-out evaluation. No audited eligible quartet or calibration of component-selective unit anchors, the Gaussian latent law or known nonlinear warp has been established. Existing construct/isoform, assay and chemical-provenance caveats remain unresolved.

An interaction contrast obtained by querying a frozen GNN describes that predictor. It does not measure a biological interaction or validate these generative assumptions. A GNN is not required for the present method. Synthetic fixtures validate neither siRNA biology nor the biological suitability of the stipulated source experiment.

**Preservation and resource audit.**

Read-only hash checks found no mismatches in the frozen learner-code commitments, existing artifact manifest, 55 protected historical/user files, or 24 prediction commitments. The linked historical ledger hash also matches. The finalized learned-target ledger SHA-256 is `29e33d348dbf42e18d3cf4dbe3ae34496af584f523ef5e868b0b87578a8dcc4a` and was not appended to during this review.

The completed pilot recorded 814.115899 CPU seconds, 817.592598 aggregate job-wall seconds and 970.392598 conservative charged seconds; zero GPU use; one numerical worker, thread and core affinity. The saved affinity-violation records are empty. The intentional timeout is retained, as is the disclosed initial historical-ledger append of .101903 charged seconds. That append prevents claiming literal zero historical-ledger writes during the original pilot, despite otherwise successful preservation checks. No remaining allowance was treated as authorization for this review to train. These are the original pilot's accounting values, not a claim that reading and writing this review consumed zero CPU.

Evidence entry points: frozen configuration (omitted repository reference: `../../../runs/interaction_recoverability/learned-target-v1-20260912-164000/frozen_config.json`), raw results (omitted repository reference: `../../../runs/interaction_recoverability/learned-target-v1-20260912-164000/evaluation/raw_results.csv`), diagnostics (omitted repository reference: `../../../runs/interaction_recoverability/learned-target-v1-20260912-164000/evaluation/diagnostics.json`), [figures and prior summary](results.md), [source matrix](novelty.md), saved correctness report (omitted repository reference: `validation.md`), final resource snapshot (omitted repository reference: `../../../runs/interaction_recoverability/learned-target-v1-20260912-164000/resource_snapshot-28.json`), and unchanged ledger (omitted repository reference: `../../../runs/interaction_recoverability/learned-target-v1-20260912-164000/ledger.jsonl`).
