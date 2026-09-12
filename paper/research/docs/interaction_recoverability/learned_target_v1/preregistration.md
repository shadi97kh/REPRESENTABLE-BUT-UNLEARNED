> Publication copy of `docs/interaction_recoverability/learned_target_v1/preregistration.md`. Mathematical and result content is retained; execution material is omitted. This dated record retains its original conclusions and limitations, including conclusions superseded by later stages. Original SHA-256: `50dd9f3cfed09de1de550dd2b680c649c834e240c752da5ddab0d92cd2dcf52b`.

# Frozen pilot and measurement rules

Status: frozen and executed; numerical outcomes belong in `results.md`.

The authoritative machine-readable configuration is [the frozen config](../../../../../configs/interaction_recoverability_learned_target_v1.json). Its SHA-256 is `567e577160f03f38ad7c2b4b8b604cccfa73efcd58b9e103dd47be554b67b1e6`. All learned-target Python source hashes are committed there. The freeze occurred before final source generation and before final target access.

One complete development dataset took 30.99345 CPU seconds / 34.73277 wall seconds inside the fitter. Its full-pilot projection was `24*1.5*34.73277 = 1250.37961` seconds, within the fresh allowance after reserving 600 seconds and an additional 60 seconds. The **full** plan was therefore chosen by the predeclared runtime rule, without using development target errors for that choice.

| Dimension | Frozen choice |
|---|---|
| Independent synthetic families | Oscillatory and localized analytic perturbations of exp(x), both globally admissible |
| Separation delta | .1 and .03 |
| Sample size n | 256 and 1024 **per source stratum**; total sample size 5n |
| Replicates | Three independent datasets per family/delta/n cell; 24 final datasets |
| Truth shifts | Four independent Uniform(.6,.9) shifts per dataset, private evaluator seed |
| Neural family | Two width-four tanh modulation profiles, jointly fitted four shifts |
| Source optimizer | Adam .03, checkpoints 100/200/300, CRPS64; no target-driven restarts |
| Source selection | Minimum C-role CRPS checkpoint; pooled controls use step 300 |
| Correction selection | Three lambdas; explicit heuristic score-based proxy, no E or target tuning |
| Enrichment | At most four rounds, tau .001, independent 36-direction analytic search bank |
| Repeated roles | Three cyclic, stratified, disjoint F/C/E parts |
| Direct-target control | Same observations/class information, constrained basis sieve, 40 steps per endpoint |
| Primary loss | Ordinary squared error of the original-scale interaction; no clipping or outcome transform |
| Screening threshold | At least 10% aggregate MSE reduction versus strongest **complete** non-oracle control |

Ten prediction methods are frozen: neural shared plug-in, neural fixed/adaptive/random/rich corrections, basis shared plug-in/rich correction, neural pooled plug-in, basis pooled plug-in, and basis direct-set midpoint. Null/failed predictions remain raw rows and are excluded from complete-control eligibility. The latter exclusion cannot justify superiority to an uncomputed optimal control. All optimization traces and numerical failures remain saved.

Three replicates per cell support descriptive pilot comparisons, **not** precise confidence intervals, rate fitting, reliable significance tests or a benchmark campaign. The strongest complete control is selected for the final descriptive comparison after predictions are frozen; it never feeds back to training. Pairwise errors, sample counts, all failures, and both total and fully charged per-method CPU costs are retained.

Stops are the compute reserve, a nonfinite/failed fit (with existing artifacts retained), and completion of the frozen dataset list. There is no adaptive expansion to additional seeds, separation values, methods or tuning to obtain a win. Source-law feasibility failure in the direct-set control is an outcome, not an instruction to restart. A theorem gap is retained even if a finite-sieve method wins.

The synthetic truths are defined in the frozen evaluation-only module. Their modulations obey the original derivative envelopes, but neither family spans the whole analytic class. Joint target integrals use independent 256/512-node refinements on `[-12,12]`. Evaluation-only targets and family labels are not passed to the learners. All configurations, learner outputs and per-fold roles can be reconstructed from the saved configurations, learner outputs, and run artifacts.
