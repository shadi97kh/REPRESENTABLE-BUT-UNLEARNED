"""Export supported historical tables and the P4 decision without fitting."""
from __future__ import annotations
import csv
import datetime
import io
import json
from pathlib import Path
from .p4_resources import checked_record, sha, usage, usage_delta, write_record


def latest(pattern):
    paths = sorted(p for p in Path('runs/revision').glob(pattern) if not p.name.endswith('-rows.json'))
    if not paths:
        raise FileNotFoundError(pattern)
    return paths[-1]


def write_csv(path, rows):
    with Path(path).open('x', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with Path(str(path)+'.sha256').open('x') as f:
        f.write(f'{sha(path)}  {Path(path).name}\n')


def main():
    start = usage()
    paths = {k:latest(v) for k,v in dict(audit='p4_audit-*.json',dryrun='p4_dryrun-*.json',
                                        cycle='p4_cycle-*.json',tests='p4_tests-*.json').items()}
    recs = {k:checked_record(p) for k,p in paths.items()}
    if recs['tests']['exit_code'] != 0:
        raise RuntimeError('latest contract tests did not pass')
    audit, smoke = recs['audit'], recs['dryrun']['smoke']
    data = checked_record(audit['row_record'])
    h, d = audit['historical'], data['summary']
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%d-%H%M%S-%f')
    prefix = Path('runs/revision')/f'p4_export-{stamp}'
    files = {}
    # No invented row-level prediction, paired group CI or measured selection.
    p3 = checked_record(audit['config']['historical_records']['p3_pilot'])
    prediction_rows = []
    for alias in h['p3_aliases']:
        arm = alias['names'][0]
        for seed, rho in enumerate(p3['summary'][arm]['spearman_per_seed']):
            prediction_rows.append({'evidence':'historical_exploratory_summary_ridge',
                                    'arm':arm,'actual_class':alias['actual_class'],
                                    'features':alias['feature_dim'],'seed_changes_split':seed,
                                    'spearman':rho,
                                    'paired_delta':h['p3_full_precision_contrasts'][arm]['paired_deltas'][seed]})
    path = str(prefix)+'-historical_prediction.csv'
    write_csv(path,prediction_rows)
    files['historical_prediction'] = path
    p2 = checked_record(audit['config']['historical_records']['p2'])
    completion = []
    for cap in [.005,.01,.02,.05,.1,.2,.5,1.,5.]:
        for arm in p2['aggregates']['overall']['arms']:
            if arm == 'lp_polytope_external':
                continue  # its degree constraint is invalid for layered paths
            for rung in range(4):
                indexed = [(c,r) for c,r in zip(p2['cases'],p2['raw']) if c['case']['rung']==rung]
                n = sum(c['arms'][arm]['decision_certified'] and r[arm]['profile']['wall_s'] <= cap for c,r in indexed)
                completion.append({'evidence':'historical_P2_completed_updates_only_no_learned_arm',
                                   'arm':arm,'rung':rung,'budget_s':cap,'n':len(indexed),
                                   'completed_decisions':n,'completion_only_coverage':n/len(indexed)})
    path = str(prefix)+'-historical_coverage.csv'
    write_csv(path,completion)
    files['historical_coverage'] = path
    path = str(prefix)+'-target_counts.csv'
    write_csv(path,[{k:v for k,v in r.items() if not isinstance(v,list)} for r in data['per_gene']])
    files['target_counts'] = path
    used = []
    for p in sorted(Path('runs/revision').glob('p4_*.json')):
        r = checked_record(p)
        if 'resources' in r:
            used.append({'path':str(p),**r['resources']})
    cpu = sum(r['process_cpu_s']+r['child_cpu_s'] for r in used)
    wall = sum(r['wall_s'] for r in used)
    ledger = audit['resource_ledger']
    table = [
        ['Actual architecture / nonadditivity','PASS, no-fit diagnostic',f"CertifiedModel; {smoke['distinct_valid_adjacencies']} valid graphs; max additive-fit residual {smoke['fixed_model_nonadditivity']['max_abs']:.6g}; exact gamma=0; residual gradients nonzero."],
        ['Independent nonlinear capacity fit','BLOCKED',f"Identifiable control implemented; equal-marginal label contrast >= {smoke['equal_marginal_label_difference_min']:.3f}. Zero optimizer steps executed."],
        ['Predictive noninferiority','NOT EVALUATED','Predeclared Spearman margin 0.02; candidate activity/NDCG loss <=0.02. P3 was summary ridge only.'],
        ['A+B / learned-bound useful gain','NOT EVALUATED','P1 repaired A-init but pairwise ridge won; P2 did not train a policy.'],
        ['Learned versus fixed, equal total wall time','NOT EVALUATED','New comparison charges setup, queries, policy inference and offline amortization; training paths not end-to-end validated.'],
        ['RNA decisions versus measured activity','LIMITED','Historical pooled ranking exists; no stored predictions for selected-activity or paired group resampling.'],
        ['Numerics / RNA context','CONDITIONAL',f"float_unverified. No verified construct; {d['corrected_conditional_intersection_with_p3']} retained historical rows under a stated shared-insert assumption."],
        ['Defensible operator novelty','UNESTABLISHED','Exact family support composed with known relaxations and learned query allocation; no first/new-mechanism claim.'],
        ['Resources / remaining','UNVERIFIED BALANCE',f"Instrumented audit/check use before export: {cpu:.3f} CPU seconds, {wall:.3f} wall seconds, zero GPU. Actual prior balance unknown; spendable verified training allowance is zero."],
    ]
    note = ['# P4 — actual-model capability and useful-certification decision', '',
            '**Decision: P4 is not scored because remaining training resources cannot be verified.** The attachment explicitly permits the read-only/implementation fallback. GO, ML GO / RNA LIMITED, and NO-GO require experimental outcomes; an unexecuted learned comparison is neither a success nor a failure. No training was launched.', '',
            '| Question | Decision | Evidence |','|---|---|---|']
    note += ['| '+' | '.join(row)+' |' for row in table]
    note += ['', '## What changed the evidence', '',
             'P1 evidence is the hash-checked **20260905-223008** revision, not the older 0.7253/0.4573 pilot. A-initialized fine-tuning improved size-holdout MSE from 10.6973 to 9.9644 on the complete-matching interaction task; pairwise ridge was effectively exact there. On layered-DAG interactions, pairwise ridge had MSE 1.7952 versus 2.6356 for the residual. These are per-structure synthetic labels, so they do not answer the new ensemble-level capacity question. The P1 function named `path_matching` actually enumerates matchings on the complete graph; P2 uses adjacent-edge path matchings.', '',
             '| Stage | Actual architecture / structural dependence | Baselines and labels | Budget comparison |',
             '|---|---|---|---|',
             '| P1 S4–S5 | Actual CertifiedModel; fitted residual nonadditive but weak | Additive/pairwise ridge and neural controls; synthetic per-member labels | Paired seeds and size holdout; no learned-bound comparison |',
             '| P2 | Random initialized and perturbed CertifiedModel per case, not a trained P1 checkpoint | Box/combined, conditional/selective, LP, exhaustive reference; no efficacy labels | Seven arms, 108 cases; old timer defect below |',
             '| P3 | Three ridge functions behind six names; sampled edge count is additive | Guide/site onehot and four summaries; unclipped reporter labels | Five different inspected partitions; no actual GNN or verifier |',
             '| P4 | Actual-model forward/backward and ensemble adapter checked | Independent higher-order control defined, not fitted | Read-only/no-fit work only; automatic run refusal |', '',
             '**P2 timing correction.** The live timer is read again after other arms finish. Immediate durations already stored in the record give these corrected medians; the historical results are untouched.', '',
             '| Arm | Reported wall s | Own immediate wall s |','|---|---:|---:|']
    for arm,v in h['p2_timing_reconciliation'].items():
        note.append(f"| {arm} | {v['median_reported_later_wall_s']:.6f} | {v['median_immediate_own_wall_s']:.6f} |")
    note += ['', 'P2 rung 0 remains **16/27 = 59.26%** coverage for box, fixed conditioning and deterministic selective. Across all rungs, coverage is 46/108, 54/108 and 53/108 respectively. Deterministic selective is faster than always conditioning but slower than box by its own immediate timer. The old claim that it was faster than box is withdrawn. Its learned-policy gate must be reassessed under a correct matched-wall comparison. The LP degree≤1 constraint excludes degree-2 internal vertices of layered-DAG paths, so that arm is not a valid general outer relaxation there. No numerical soundness proof follows from the recorded zero violations.', '',
             '## P3 claims-to-evidence corrections', '',
             '| Historical claim | Corrected evidence |','|---|---|',
             '| Six additional eligible arms | Two to six means **four** additional nominal arms. |',
             '| Six architecture comparisons | Three fitted functions: 85-coefficient guide/site ridge; 89-coefficient summary ridge; 90-coefficient sampled-count summary ridge. No chemistry features and no MFE/GNN training. |',
             '| Exact identical predictions | Only metric vectors were saved. Their hashes match for aliases; prediction hashes cannot be recovered without an explicitly charged reconstruction. The cached-count aliases depend on configured execution order. |',
             '| Gene-level hull proves native insert context | UBE2B, UBE2C, UBE2E3 and UBE2V1 mix first-hit transcript coordinate systems. New hulls are per accession and transcript, conditional on one shared native insert; dense tiling does not verify that assumption. |',
             f"| 1,156 interior / 1,153 pilot | Three UBE2E3 windows are truncated by transcript ends; rows listed below. Full accounting: {d['full_rows']} unique guides/rows, {d['legacy_located_rows']} located rows over 18 targets, 1,153 pilot rows over 18 targets. The five-gene prose table was illustrative, not the complete population. |",
             f"| 906 isoform-ambiguous rows | On the old selected coordinates, 700 have identical local windows and 206 differ. Under consistent 19-mer coordinates: 1,013 historical rows have one distinct admissible local window, 117 have multiple, and 23 have none. Intersection retained: {d['corrected_conditional_intersection_with_p3']}. |",
             '| Clean novel-sequence split | No exact duplicates; substantial shared 13-mers and overlapping sites across all five partitions. Zero gene×minimum-13-mer identifier leakage was insufficient. |',
             '| Structure effect / seed-noise ratio | Withdrawn. Seeds change splits, not independent biology. Ridge has no initialization; the seed-0 graph-count cache is reused. Full-precision paired score deltas are preserved below. |',
             '| Prior sensitivity versus arm performance | Incommensurate comparison withdrawn. Old sensitivity uses theta=0.5×reported temperature, while fitted summaries use theta=0. SE is in base-pair counts, not predicted activity. It is not prediction uncertainty. |', '',
             'The corrected isoform rule maximizes the number of distinct guides mapping to one transcript within an accession, breaking ties by transcript ID; repeated positions break by coordinate. It never sees efficacy labels. Alternative-window and identical-window-only analyses are implemented but **not run**. No uniform transcript mixture is called biological expression. All contexts remain `inferred_contiguous_native_insert` or `unresolved`; `verified_construct` has zero rows. F9 stays quarantined. Primary assay identity remains YFP/H1299/48 h. Labels above one are preserved (134; maximum 1.341). Raw assay-replicate identifiers are absent.', '',
             '| Excluded row | Transcript | Start (0-based) | Transcript nt | Actual window nt |','|---|---|---:|---:|---:|']
    for r in data['three_short_exclusions']:
        note.append(f"| {r['row_id']} | {r['legacy_transcript_id']} | {r['legacy_pos']} | {r['legacy_transcript_len']} | {r['actual_window_nt']} |")
    note += ['', '| Seed | Test rows | Test guides sharing ≥13 contiguous bases with train | Cross-gene shared-13-mer pairs | Test rows with overlapping mapped sites |','|---|---:|---:|---:|---:|']
    for r in data['legacy_cross_partition']:
        note.append(f"| {r['seed']} | {r['n_test']} | {r['shared_13mer_test_rows']} | {r['cross_gene_shared_13mer_pairs']} | {r['overlap_test_rows']} |")
    for arm in ['mfe_gnn','sampled_ensemble_gnn']:
        v = h['p3_full_precision_contrasts'][arm]
        note += ['', f"`{arm}` minus guide ridge: mean **{v['mean_delta']:.16f}**; paired deltas `{v['paired_deltas']}`. These are split-paired descriptive differences, not a group-bootstrap confidence interval."]
    note += ['', '## Implemented, run, and bounded blocker', '',
             f"The independent control uses a prespecified two-of-three motif interaction on an enumerated 26-member family, with additive and pairwise ridge controls. Paired distributions have equal first-order marginals to {smoke['equal_marginal_pair_max_error']:.3g} but expected-label differences ≥0.16. The target is not exactly pairwise (projection residual max 0.2). These are algebraic checks, not a capacity-fit pass. Actual-model nonadditivity max is {smoke['fixed_model_nonadditivity']['max_abs']:.6g}. Ensemble prediction equals the explicit weighted mean; loss-of-mean and mean-of-loss differ by {smoke['difference_is_structure_variance']:.6f}, demonstrating why they must not be interchanged.", '',
             f"Measured dry-run estimate for the **single** 200-anchor/400-residual-epoch control is {smoke['estimated_control_training_wall_s']:.1f} seconds on one CPU thread (four times the measured full-family forward/backward time; excludes verifier teacher generation). Any executed cycle would also be bounded by the smaller of 900 seconds and half the verified remaining CPU/wall/instance allowances. GPU use is disabled. The existing approval is 8 CPU core-hours, 1 GPU-hour, 2 wall-hours and 600 instances. Just both original pilots plus P3 could have consumed {ledger['counterexample_16_core_hours_pilots_plus_p3']:.3f} core-hours if all 16 CPUs were active: the saved wall proxies cannot prove a positive balance. This is not a claim that the budget was actually exhausted. The user confirmed no reconciled ledger is available.", '',
             'Implemented code includes the actual ensemble training objective, ridge controls, global sequence grouping, RNA model adapters, target-level ranking and paired resampling, a frozen-predictor selective policy comparison, total-time curves, and the refusal/timer guards. The **training paths remain unexecuted and are not end-to-end validated**. A full RNA run is also checked against the old instance cap; it cannot silently treat 1,130 rows as fitting within 600 instances. No fitted prior/isoform prediction sensitivity, gradient-sampling variance result, learned-verifier gain, or full RNA comparison is claimed.', '',
             'The three requested main-figure datasets are therefore not all supported. Exported historical prediction and completion-only P2 curves are marked historical/supplemental; the latter uses only completed recorded updates and is not a learned-versus-fixed benchmark. Per-target counts are available, but a measured selected-activity figure is withheld because historical predictions were not saved.', '',
             '## Prior work and Davis', '',
             'The implementable distinction is **the constrained support operator** `h_F(a + gamma*c)` over a declared graph family, combined with an exactly bypassable additive anchor and a learned choice of which conditional-support computations to request. This is a composition of existing ideas, not an established new tightening mechanism. Hojny et al. already condition topology bounds on fixed edges (§3.4, Eqs. 21–23); their SCIP-MPNN is the closest verified operator-level reference. General neural branching is prior work in Lu–Kumar. Non-crossing-family *membership* is easy to check; the nontrivial DP here is constrained optimization/counting, so the old novelty note’s “hard membership” rationale is withdrawn. No “first” claim is earned. [Hojny et al.](https://arxiv.org/html/2402.13937v2), [SCIP-MPNN code](https://github.com/christopherhojny/SCIP-MPNN), [Lu–Kumar](https://arxiv.org/abs/1912.01329).', '',
             'Davis remains absent. The required file is `gkaf479_supplemental_files.zip`, containing Supplementary Table S1, to be placed at `data/davis2025/gkaf479_supplemental_files.zip`. The normal publisher issue page → article → Supplementary data route was inspected; it redirected to an unusable minimal article page. PMC presented reCAPTCHA; it was not bypassed or parsed. Targeted UMass institutional, GitHub and Zenodo searches found no verified supplement copy. No authors were contacted. The archive was not acquired, so no S1 table version or covariate audit is claimed. Chemistry transfer, native/reporter cohorts, doses/times and replicates remain unresolved; GSE231101 is not substituted. This does not block the independent ML control. [Publisher route](https://academic.oup.com/nar/issue/53/12), [article](https://academic.oup.com/nar/article/53/12/gkaf479/8171869).', '',
             '## Reproducible commands and records', '', '```bash',
             'python -m scmp.revision.p4_decision --help',
             'python -m scmp.revision.p4_decision --audit',
             'python -m scmp.revision.p4_decision --dry-run',
             'python -m scmp.revision.p4_checks',
             'python -m scmp.revision.p4_decision --run  # expected exit 2: unverified resources',
             'python -m scmp.revision.p4_export', '```', '',
             'Each invocation writes new timestamped outputs and SHA256 sidecars. Existing model code, readout scale, checkpoints and historical records are preserved. The first test attempt had one test-harness assertion failure (it expected no outer timer while the test runner correctly maintained one); that record is retained. The corrected run passed all 54 tests with no optimizer steps.', '']
    note += [f"- {k}: `{p}`" for k,p in paths.items()]
    note += ['', '**Single recommended next action:** establish a verifiable, explicitly scoped allowance for the one bounded capability control. Do not scale RNA or verifier experiments from the historical summary null.']
    report = Path('docs')/f'p4_actual_model_decision-{stamp}.md'
    report_text = '\n'.join(note)+'\n'
    with report.open('x') as f:
        f.write(report_text)
    with Path(str(report)+'.sha256').open('x') as f:
        f.write(f'{sha(report)}  {report.name}\n')
    out = {'kind':'p4_export','report':str(report),'report_sha256':sha(report),
           'input_records':{k:{'path':str(p),'sha256':sha(p)} for k,p in paths.items()},
           'datasets':{k:{'path':p,'sha256':sha(p)} for k,p in files.items()},
           'decision':'NOT_EVALUATED_RESOURCE_BLOCKED','decision_table':table,
           'resources_before_export':used,'recorded_cpu_seconds_before_export':cpu,
           'resources':usage_delta(start),'training_executed':False,
           'main_figure_status':{'prediction':'historical_summary_only','verification':'historical_completion_only_no_learned_arm',
                                 'rna_measured_selection':'unavailable_predictions_not_stored'},
           'next_action':'verify a scoped resource allowance for the single bounded control'}
    print(write_record(str(prefix)+'.json',out))
    print(report)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
