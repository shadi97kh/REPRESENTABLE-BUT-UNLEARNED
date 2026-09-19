# Codex command: focused follow-up analyses and final manuscript revision

Work in `/home/shadi/iclr2027`.

Read this instruction completely, applicable AGENTS.md, git status, and the
current audit manifest/report before editing. Execute the unblocked work,
revise the manuscript from the results, and verify the deliverables. Do not
stop at a plan, a literature assessment, or prepared code.

This instruction supersedes earlier versions of this SAME instruction file.
Continue from the completed 50-page manuscript, not the superseded 80-page
version. The prior instruction's ban on additional diagnostic analyses is
superseded only by the specific analyses authorized below. The prohibitions
on new training fits, new scripts, and unnecessary duplicate files remain.

## 1. Scope, preservation and canonical paths

ZERO new training fits. ZERO training optimizer updates.
ZERO additional first-party executable scripts, wrappers or notebooks.
No architecture search, broader hyperparameter/seed sweep, synthetic-label
training, contribution gates, polytope revival, or reachability-project work.

Use and extend the existing canonical runner:

    analysis/structural_attribution_audit/audit.py

Reuse these canonical outputs and working source:

    analysis/structural_attribution_audit/audit_results.sqlite
    analysis/structural_attribution_audit/audit_manifest.json
    analysis/structural_attribution_audit/audit_report.md
    papers/interaction_recoverability_iclr2027/current/source/
    papers/interaction_recoverability_iclr2027/current/artifacts/main.pdf
    papers/interaction_recoverability_iclr2027/current/artifacts/manuscript_source.zip

Inventory the actual current state. Preserve user edits and historical
results. The latest reported inventory is 26,601 persistent files and 743
scripts, but measure the actual values instead of assuming they still hold.
The historical fit ledger is 3,351; change that count only if a real accounting
error is demonstrated, not because old checkpoints are evaluated again.

Inspect existing functions, stages, tables, prediction exports, masks and
checkpoint metadata before implementation. Extend existing stages/functions.
If a separate stage is needed for dependency clarity, implement it inside
this SAME runner; do not introduce another executable. Keep the existing
`all` and `verify-resume` interface. Do not invent or assume other CLI flags.

Store new results as rows/tables in the same SQLite database and summarize
them in the same report. Update the same figures, manuscript, PDF and source
ZIP. No timestamp directories, version chains, new review packets, duplicate
reports, or extra source/data ZIPs. Read existing evidence in place. Reuse one
owned temporary/build area. Use PYTHONDONTWRITEBYTECODE=1.

Keep real scientific dependency hashes and cache invalidation accurate.
Changed calculations invalidate their actual dependents; editorial changes
must not trigger unrelated scientific computations. Never hide a relevant
code/input change to force a cache hit. Preserve failed executions and explain
which stages executed, which were reused and which were blocked.

If the downloaded instruction has a '(1)' or '(2)' suffix, reconcile it to
the canonical registered instruction path, preserving its superseded hash
and unrelated edits. Do not maintain two active instruction files.

No submission, public posting, git push, external contact, paid resources or
wet-lab work is requested. Public documentation/asset retrieval and the
bounded official inference/preprocessing attempt in section 6 are in scope.

## 2. Verify and integrate endpoint and marginal-effect diagnostics on B3

Question: does B3 demonstrate an interaction-specific failure of an otherwise
strong activity predictor, or broader predictive weakness on the same panel?
Use the existing source-held-out B3 prediction export and checkpoints. Do NOT
substitute the source-included deployment checkpoints used in the separate
parameter perturbation experiment.

The review reconstructed these results from the historical export:

    runs/sirna_gnn_empirical/v3-20260914T013314Z/b3_primary/interaction_predictions.csv

The copy examined in the review had SHA256:

    ebebee69d3a9db014c031e0a17b0dc9e75cd4c2157de179b191233146c9a4b55

Locate the authoritative local input through the manifest. If its identity
or numbers differ, investigate and report the difference; do not force a
match to the review or silently substitute a historical target.

Required reconstruction:

1. Identify the 165 unique measured AS/SS endpoint states by their molecular
   and assay identities. Deduplicate repeated appearances of 00, 10 and 01
   across rectangles. Do not count the 560 rectangle-corner appearances as
   560 distinct observations. Verify repeated labels/predictions agree.
2. Use the ten actual final seeds 1103, 2207, 3301, 4409, 5519, 6637, 7753,
   8867, 9973 and 11027 for the four current methods. Exclude saved 'ensemble'
   rows and deterministic-method seed-0 rows from constituent averaging.
3. Reconstruct the four ensembles. Verify their 140 interaction errors agree
   with current Table 6 on the same labels and identities.
4. Report original-scale endpoint MSE, R-squared, correlation, measured and
   predicted SD, bias, and the explicitly retrospective endpoint test-mean
   constant. Recover a training-fitted constant only from genuine matching
   training metadata; do not manufacture it from these evaluation labels.
5. Report the 14 antisense replacement effects f(A,0)-f(0,0), the 10 sense
   effects f(0,B)-f(0,0), and the original 140 B3 interactions separately.
   For each, report MSE, correlation, mean and SD of predictions/labels and
   the zero-effect comparison. Keep directions and response scale fixed.
6. Keep seed variability separate from ensemble performance. No independent-
   rectangle or independent-marginal confidence interval is justified.

Reference values from the review, to verify rather than hard-code:

| Method | Endpoint R2 | AS-effect correlation | B3 MSE |
|---|---:|---:|---:|
| Corrected GNN | 0.0318440474 | 0.0287789914 | 0.0683350733 |
| Corrected no-message | 0.0318098253 | 0.0197534598 | 0.0682810009 |
| Chemistry tree | 0.0168725623 | 0.1167793546 | 0.0719387701 |
| Token CNN | 0.032767 approximately | 0.055427 approximately | 0.068408692 approximately |

The endpoint retrospective constant MSE is 0.0798942224. For the GNN,
endpoint prediction SD is 0.0402647535 versus measured SD 0.2826556604;
AS-effect predicted SD is 0.0243522030 versus measured SD 0.1706159407.

Interpret the verified result honestly. Weak endpoints and weak marginal
effects mean this panel does not cleanly isolate an interaction-specific
failure of a strong activity predictor. Do not infer the cause of attenuation
from that observation. Do not replace the original reference-based B3 target
with a two-way centered interaction and silently treat them as equivalent.

## 3. Run one label-independent preprocessing/rank ablation

Question: which representation or preprocessing operation causes distinct
reported chemical states to merge, and which contrast restrictions follow?

The current audited representation has 165 states, 90 classes, rank(C)=140,
rank(CH)=72, 68 independent homogeneous constraints on the predicted contrast
vector, and 0/140 individually forced-zero rows. The forty audited instances
share a partition; they are not forty independent representation designs.

Hold fixed the raw panel, endpoint identities, response convention and C.
Inspect the actual pipeline and choose the following label-independent
variants before examining their projected-label residuals:

A. Canonical parsed chemical records, retaining the verified positional
   chemistry, bases, strand identities, terminal/linkage information and
   assay context that are genuinely available.
B. Actual feature tensors before training-support masks/feature selection.
C. Actual tensors after each relevant information-removing operation, with
   its saved training-only mask and method-specific consumed branches.
D. A diagnostic representation retaining the verified chemical distinctions
   removed by masking, using the existing schema or a clearly specified
   canonical representation of reported chemistry.

Do not add outcome-bearing identifiers or arbitrary row IDs as model
features. Do not invent missing geometry, linkage directions, stereochemistry
or chemical aliases. Do not pick variants based on their error reduction.
Pure raw-state identity, if used as an algebraic upper control, must be labeled
as such rather than advertised as a useful learned chemical representation.

For each distinct applicable preprocessing instance and variant, compute:

- raw/admitted state counts and invalid/unsupported counts;
- exact equivalence classes and their membership/multiplicity;
- rank(C), exact rank(CH), left-null dimension, individually zero rows;
- measured-vector projection residual and the same explicit denominator;
- which original classes split or merge and which operation accounts for it.

Reuse exact rational rank calculation and existing numerical projection
routines with declared tolerances. Deduplicate identical partitions for
computation but preserve the instance-to-partition mapping. Check the actual
ensemble partition rule instead of carrying over a constituent result when
partitions differ. Keep literal tensor equality separate from conservative
full-input checks or real-arithmetic routing equivalence.

Give two or three concrete pairs of distinct reported chemical states that
collapse, identifying the precise parsed field/tensor coordinate removed.
Record unsupported vocabulary separately from a valid encoding collision.
Check that a claimed representation refinement really refines the partition;
where refinement holds, its attainable contrast span cannot shrink.

Report the result even if the rank stays 72 or the suspected mask is not the
cause. All four current models share preprocessing, so a common restriction
cannot automatically be attributed to message passing or to published models.

This is an encoding-only ablation. Do NOT score a frozen checkpoint under
changed preprocessing as though it were a properly trained alternative.
Higher rank establishes additional representable contrasts, not better
prediction, identified biological effects or a new architecture result.

## 4. Complete the weighting/baseline comparison and bounded influence check

The manuscript's statement that the source-excluded tree beats both training
constants needs its row-weight qualifier. Current Table 22 reports:

| Evaluation weighting | Tree MSE | Training-row constant | Training-group constant |
|---|---:|---:|---:|
| Rows | 0.052931 | 0.057080 | 0.054892 |
| Equal study component | 0.060055 | 0.061373 | 0.059520 |

Thus under equal study weighting the tree is worse than the training-group
constant by about 0.000535 MSE. Recompute from full-precision matching saved
predictions/constants, not the rounded table.

Use the existing current primary-assay and source-excluded records. Report
GNN and tree minus each training-fitted constant under BOTH evaluation
weightings. Keep the fit partition, evaluated rows, labels, response convention
and weights explicit. Never replace a training constant with a test mean.

Reuse matching conditional intervals if already stored. For missing ones,
extend the existing joint paired component-resampling implementation only as
needed. For grouped scientific claims resample whole study-linked components;
keep sequence-conditional analyses separately labeled. Use the existing fixed
analysis seed and 10,000-draw convention, with prediction/selection fixed.
Preserve row-mass ratios for row weighting and equal study mass for the other
estimand. Do not treat model and baseline errors as independent samples, or
subtract interval endpoints. Report the small number of components and the
conditional nature of these intervals. No new p-value campaign or post-hoc
equivalence threshold is needed.

Perform one descriptive leave-one-study-component-out RESCORING analysis for
these comparisons. Leave the fitted models and their original training-
partition constants unchanged; omit one evaluation component and renormalize
the specified evaluation weights. Do not retrain. Report omitted component,
remaining coverage, loss differences and whether the descriptive ordering
changes. Recompute any retrospective R2 denominator on the scored subset;
flag zero/near-zero denominators. These omissions are influence diagnostics,
not new independent studies or uncertainty about a refitted learner.

Present one compact matched comparison, with the existing uncertainty that
actually supports it. Do not use ten-seed SD as between-study uncertainty.

## 5. Quantify the practical size of the existing perturbation result

Reuse the proven training-independent blocks, three saved seeds, original
and corrected arms, fixed Rademacher direction and scales 0 and +/-0.25.
No new parameter directions, larger scales, optimization, checkpoint search
or training. Do not save perturbed checkpoints.

Recheck the original proof's applicability and exact restoration/invariance
checks when the affected inference is executed. Preserve all training outputs,
zero-scale controls, corrected-gate null controls and B3 null outcomes.
Training-prediction independence is not automatically independence of model
selection, validation, regularization or every training-procedure choice.

First reuse stored full prediction/gradient vectors when available. If only
summary records were saved, run the bounded missing inference/derivatives
through the same existing runner and record that work accurately.

Add these summaries:

1. Within each of the 17 exact APP assay pools, compare baseline and perturbed
   top-five selection masses. Use the existing fractional boundary-tie rule.
   With selection masses s_i and t_i, each summing to k=min(5,n), report
   overlap sum_i min(s_i,t_i)/k. Report changed-pool counts, pool sizes and
   per-case distributions. The pools overlap and are not independent trials.
2. Compute average-rank displacement within each pool, normalized by n-1,
   with median, upper quantile and maximum. Keep global APP rank-change
   counts as a separate diagnostic rather than a candidate-selection claim.
3. Report APP absolute prediction-change distributions, RMS change, and
   change relative to the original within-pool prediction spread. Flag
   zero/near-zero reference spread instead of silently dividing by epsilon.
4. Report absolute L2 and maximum-coordinate raw-gradient changes, original
   gradient norms, and relative L2 changes where the reference norm is
   meaningfully nonzero. State the threshold before calculation and count
   undefined/near-zero cases. Use a fixed coordinate domain and scale.
5. Retain the original 16 lexicographic queries for exact comparison. If
   coverage beyond them is claimed, choose at most 128 total unique queries
   by a deterministic, response-independent scheme covering assay pools and
   training-support patterns, including those 16. Record the sampling rule
   and selected IDs before inspecting new gradient outcomes. Do not choose
   the queries with the largest observed change. Report actual coverage and
   do not present sampled maxima as maxima over all APP records.

Keep the raw-gradient definition explicit: differentiation with respect to
encoded X while adjacency, masks and context are held fixed. These are local
Euclidean sensitivities, not feasible chemical interventions. Do not infer
feature-importance changes from changes in prediction ranks, or biological
explanation validity from raw-gradient variation.

Retain null or practically small results. If selections barely change, say
so and reduce the abstract's emphasis accordingly. The prior APP maximum
0.003300 is approximately 0.33 percentage points on the fractional scale;
large global rank-change counts alone do not establish meaningful impact.

## 6. Make one bounded attempt through the official ENsiRNA route

This is the only new external implementation route requested. Do not restart
an open-ended baseline search or repeat known failed MEG-mod/ModMapper paths
without materially new official assets.

The official repository currently documents a Docker route and ENsiRNA-mod
inference:

    https://github.com/tanwenchong/ENsiRNA

Treat the documented `tanwenchong/ensirna:v2` image and ENsiRNA-mod
`easy_run.py` as leads to VERIFY against the live official instructions.
A tag is mutable; record the actual image digest, code version and assets.
Do not assume the image resolves the previously missing Rosetta/PDB inputs.

1. Check whether this route was already attempted. Reuse any cached image,
   weights, geometry and successful environment instead of duplicating them.
2. Inspect documented requirements, available runtime, storage and GPU
   compatibility before a pull. Use only already available authorized local
   compute; do not install a privileged daemon, change host permissions,
   disrupt another job, or acquire a paid resource.
3. Limit this availability/example preflight to one official route and one
   documented example, with a 20-minute wall-time cap. Record timeout or
   missing capability as the concrete blocker. Inspect pull size and disk
   capacity first; do not begin a pull whose known requirements do not fit.
4. If the official example runs, inspect the complete consumed preprocessing
   for the fixed B3 panel and the already declared conformance pairs. Reuse
   endpoint states. Cap this additional fixed-panel audit at 20 minutes;
   report completed/unsupported/blocked cases with explicit denominators.
5. Include geometry, sequence embeddings, atom IDs, tokens and every other
   consumed branch in a full-pipeline equivalence claim. A matching chemistry
   fingerprint alone remains a partial-branch result.
6. Determine training/validation overlap before reporting predictive scores.
   A checkpoint trained on Bramsen is ineligible as held-out B3 performance
   evidence. Unknown overlap stays unknown. A complete encoding audit can
   still be reported with its proper scope without predictive claims.

ZERO published-model retraining is authorized in this follow-up. If faithful
predictive comparison requires retraining or unresolved scientific inputs,
record the exact remaining requirement and continue the other work. Do not
invent a geometry/token representation, patch scientific semantics to obtain
a score, or claim baseline reproduction from a help-command import.

Record any retained third-party assets separately from first-party scripts
and scientific outputs. Store logs and metadata through the existing report,
manifest and results store. Do not repull or rerun this on unchanged resume.
A failed/blocked external preflight does not stop sections 2-5 or the paper.

## 7. Correct the actual manuscript and align the claims with the results

Retain the current paper direction and title:

    Auditing Chemistry-Effect Claims in Modified-siRNA Prediction

Retain the central distinction between aggregate activity performance,
measured chemistry effects, input representation and training support.
Make no general XAI-validity theorem, novel architecture, biological
validation, or universal benchmark-failure claim.

Corrections identified in the DELIVERED 50-page PDF, regardless of what the
previous completion report asserted:

- Page 44, Appendix A11.7 still says: 'A predictor that had learned source
  levels would show the opposite sign there.' Remove that inference and any
  semantic equivalent. Report observed evaluation covariance only.
- Page 47, Table 23's caption generalizes negative between-source covariance
  to all methods. Qualify it to the GNN; the displayed tree values are
  positive. Check every method/scenario against the actual stored values.
- Qualify 'the tree beats both deployable constants' with 'under row
  weighting'. Present the equal-study result where it changes the ordering.
- In the abstract, replace '68 independent combinations of the measured
  contrasts are unattainable' with '68 independent homogeneous constraints
  on the predicted contrast vector'. Rank counts constraints; it does not
  establish 68 separately nonzero measured discrepancies.
- Replace unqualified 'ranks' with 'global prediction ranks' for the original
  experiment. Label any new within-pool selection results separately.
- Narrow 'source-level effects do not explain the pattern' to the recorded
  GNN covariance decomposition. Negative evaluation covariance does not
  establish absence of source-dependent learning.
- Page 50's command lost its spaces. Render a genuinely copyable command:
  `python analysis/structural_attribution_audit/audit.py all`.

Integrate the endpoint/marginal comparison and preprocessing-rank ablation
without framing weak endpoint prediction as a demonstrated failure specific
to interactions. Preserve the original reference-based B3 result and its
one-background/shared-measurement limitations. Preserve the B2 bundle's four-
position confounding and unresolved comparison with simple controls.

Preserve the verified decomposition for the ORIGINAL audited representation:

    0.06833507332376675
      = 0.009979515952211788 + 0.058355557371555

The residual is 14.6317% of recorded uncentered contrast energy and 14.6038%
of fitted GNN error; 85.3962% of the latter is in the permitted subspace.
The unrestricted decoder need not be realizable by the fitted model family.
Neither this decomposition nor new rank ablations identify the causal
mechanism of attenuation or remove measurement noise.

Preserve the original 0/140 cancellation finding, source-held-out identities,
protocol dependence, unfavorable/null outcomes, 299 qualified discrepancies,
S1 quarantine, unresolved dose/assay metadata and unknown shared-control
covariance. New published-method outcomes get their own version/scope;
do not rewrite the history of the pinned partial checks.

Rebalance the abstract toward the actual audit findings. The rebuttal of an
earlier mistaken interpretation of pooled R2 should not dominate the claimed
contribution. Support any claim that published work routinely makes a
particular inference with a concrete primary-source example, or soften that
premise. Do not count elementary rank-nullity or variance identities as new
mathematical theory.

## 8. Preserve the page, figure, source and evidence requirements

The ENTIRE final PDF must be AT MOST 50 physical pages, including statements,
references and appendix. Keep exactly nine main-text pages in the existing
verified official anonymous ICLR 2027 style. Keep the architecture/audit
workflow on page two and at least three useful results/workflow figures in
the main text, updating existing assets in place. Do not pad to 50 pages.

Make space by replacing repetition or less relevant historical summaries.
Keep raw predictions, per-seed dumps, execution logs, full witnesses and hash
lists in their existing evidence stores. Preserve all definitions,
assumptions and complete proofs required for the retained claims inside
the manuscript. Preserve historical evidence outside the current source.
Do not shrink official fonts, margins or spacing, omit unfavorable evidence,
or create another raw-data PDF to evade the cap.

The editable source and source ZIP must have exactly four top-level entries:

    main.tex
    appendix.tex
    references.bib
    figures/

Only final referenced artwork belongs in figures/. Keep styles, builds,
compiled manuscript PDFs, raw data, scripts and reports outside the source.
Keep figure and table captions below their objects; standalone artwork has
no embedded caption. Generate numerical tables from verified stored results.

Update the existing report's claim/evidence and proof-dependency map. Fix
stale appendix references and captions that no longer match compressed
contents. Search normalized EXTRACTED FINAL PDF text as well as source;
checking only a generator template did not remove the page-44 sentence.

## 9. Verification, execution accounting and delivery

Use meaningful checks in the existing runner, not new test scripts:

- exact endpoint deduplication and current-table reproduction;
- input/partition identities, complete-input equivalence and exact ranks;
- projection residuals, orthogonality and original-model in-span checks;
- correct weighted paired differences and uncertainty units;
- fractional tie selection masses, rank domains and gradient query coverage;
- training-prediction invariance, zero controls and restored checkpoint hashes;
- no training optimizer construction/updates or new fit records;
- concrete scope/overlap accounting for the official external attempt;
- unchanged-resume scientific stages and file inventory.

Generate report/receipt accounting from one explicit execution-record snapshot,
with post-snapshot report/build/finish costs identified separately. Retain
failed calls. Separate wall, process CPU, subprocess CPU, GPU inference if any,
peak memory and unmetered administration; do not double-count overlapping
measurements or label the entire task as zero compute.

Compile the final manuscript and inspect every page. Verify actual page
count/main boundaries, anonymity, cross-references, citations, equation/figure
labels, captions and readable layouts. Independently compile the extracted
four-entry ZIP with the existing external official style dependencies in one
reusable build area. Verify source identities and all-page text agreement.
Check that the delivered PDF, rather than only the source, contains the
scientific corrections above.

After implementing this specification, execute the known existing interface:

```bash
cd /home/shadi/iclr2027
export PYTHONDONTWRITEBYTECODE=1
python analysis/structural_attribution_audit/audit.py all
python analysis/structural_attribution_audit/audit.py verify-resume
```

These commands must include the newly implemented analyses through `all`.
Running the old unmodified runner alone would only reproduce the old paper.
On an unchanged rerun, completed or terminally blocked external stages must
not redownload assets or repeat scientific work. A later genuine unblocking
input can invalidate the appropriate stage through the normal dependency
mechanism. Unchanged resume must add zero persistent files and execute zero
scientific stages; report bookkeeping in its documented scope.

Return a concise completion report with:

1. B3 endpoint, marginal and interaction diagnostics, actual identities and
   whether the reconstructed review values agree.
2. Stagewise encoding class/rank/residual results and concrete collision
   examples identifying the responsible preprocessing operation.
3. Both-weighting model/constant differences, relevant conditional intervals
   and leave-one-study influence outcomes.
4. Within-pool selection/rank changes, gradient magnitudes and actual query
   coverage, including nulls and all prespecified perturbation cases.
5. Official-pipeline result: verified environment/example, full versus partial
   encoding scope, training overlap, counts and precise remaining blockers.
6. Corrected claims, retained limitations and actual PDF page boundaries.
7. Zero new fits/optimizer updates/additional scripts, before/after counts,
   resource snapshot, actual rerun stages and unchanged-resume verification.
8. Exact tested execution/resume commands and links to the SAME canonical
   PDF, source ZIP, report, SQLite store and manifest.

Complete every unblocked item. Do not stop at the first external blocker,
claim a stage was executed when it was merely cached, or propose another
training campaign when the evidence contradicts a preferred interpretation.
