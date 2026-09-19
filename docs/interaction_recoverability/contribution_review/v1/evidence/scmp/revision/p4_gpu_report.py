"""Export executed P4 evidence, without fitting or changing selection criteria."""
from __future__ import annotations
import argparse
import csv
from pathlib import Path
import numpy as np
from .p4_resources import checked_record, sha, write_record, usage, usage_delta


def fmt(x, digits=4):
    return 'unavailable' if x is None else f'{x:.{digits}f}'


def main(cycle_path):
    start = usage()
    cycle = checked_record(cycle_path)
    assert cycle['status'] == 'COMPLETE_ALL_REQUESTED_COMPARISONS'
    for item in cycle['artifacts'].values():
        assert sha(item['path']) == item['sha256']
    def read(name):
        return checked_record(cycle['artifacts'][name]['path'])
    cap, rna, ver, inputs, pre = [read(n) for n in ['capability','rna','verifier','rna_inputs','prerun']]
    cfg, rows, test = pre['config'], inputs['rows'], rna['split']['test']
    outprefix = Path(cycle_path.replace('-cycle.json',''))
    exported = {}
    def csv_out(name, data):
        path = Path(str(outprefix)+'-'+name+'.csv')
        with path.open('x',newline='') as f:
            writer = csv.DictWriter(f,fieldnames=list(data[0]))
            writer.writeheader(); writer.writerows(data)
        Path(str(path)+'.sha256').write_text(f'{sha(path)}  {path.name}\n')
        exported[name] = {'path':str(path),'sha256':sha(path),'rows':len(data)}
    def macro(comp, key, side='proposed'):
        vals = [v[side][key] for v in comp['per_target'].values() if v[side][key] is not None]
        return float(np.mean(vals)) if vals else None
    predictive = []
    for arm, comp in rna['comparisons'].items():
        m = comp['proposed']
        predictive.append({'arm':arm,'reference':rna['reference_selected_on_development'],'test_rows':m['n'],
           'spearman':m['spearman'],'delta_spearman':comp['delta_spearman'],
           'delta_ci95_low':comp['group_paired_ci95'][0],'delta_ci95_high':comp['group_paired_ci95'][1],
           'mse':m['mse'],'macro_target_ndcg_at_5':macro(comp,'ndcg_at_k'),
           'macro_target_precision_at_5':macro(comp,'precision_at_k'),
           'macro_selected_activity':macro(comp,'selected_activity'),
           'macro_selection_regret':macro(comp,'selection_regret'),
           'context_status':'inferred_contiguous_native_insert','evaluation_status':'exploratory'})
    csv_out('figure1_prediction_ranking',predictive)
    coverage = []
    for r in ver['equal_total_cost']:
        for arm, v in r['coverage'].items():
            online = next(s for s in ver['summary'] if s['arm']==arm and s['stratum']=='original_family' and s['budget_s']==r['budget_s'])
            coverage.append({'arm':arm,'stratum':'original_family','total_budget_s':r['budget_s'],
              'amortization_instances':r['amortization_instances'],'decision_coverage':v,'n':online['n'],
              'offline_charge_s':ver['offline_teacher_and_fit_wall_s']/r['amortization_instances'] if arm=='learned_selective' else 0.,
              'online_mean_gap_when_available':online['mean_gap_when_available'],
              'online_gap_unavailable_n':online['gap_unavailable_n'],'online_timeouts':online['timeouts'],
              'online_wall_p90':online['wall_quantiles'][1],'numerical_status':'float_unverified',
              'predictor_capability_gate':cap['capability_passed'],
              'scope':'diagnostic; failed capability gate precludes useful-verification claim' if not cap['capability_passed'] else 'original control family'})
    csv_out('figure2_coverage_total_compute',coverage)
    genes = sorted({rows[i]['gene'] for i in test})
    per_target, candidates = [], []
    pred = np.array(rna['predictions']['actual_A_plus_B'])
    ref = np.array(rna['predictions'][rna['reference_selected_on_development']])
    variants = {f'prior_T_{t}':np.array(r['predictions']) for t,r in rna['prior_sensitivity'].items()}
    variants['alternative_isoform'] = np.array(rna['isoform_alternative_predictions'])
    variants['independent_64_draw'] = np.array(rna['sampling']['independent_draw64_predictions'])
    variants['draw_16'] = np.array(rna['sampling']['draw16_predictions'])
    for gene in genes:
        ix = np.array([i for i in test if rows[i]['gene']==gene])
        ranks = {name:{int(i):int(rank+1) for rank,i in enumerate(ix[np.argsort(-v[ix],kind='stable')])}
                 for name,v in {'actual_A_plus_B':pred,'reference':ref,**variants}.items()}
        comp = rna['comparisons']['actual_A_plus_B']['per_target'][gene]
        per_target.append({'gene':gene,'n':len(ix),'prevalence_at_0_70':comp['proposed']['threshold_prevalence'],
             'ndcg_at_5':comp['proposed']['ndcg_at_k'],'reference_ndcg_at_5':comp['reference']['ndcg_at_k'],
             'selected_activity':comp['proposed']['selected_activity'],'reference_selected_activity':comp['reference']['selected_activity'],
             'precision_at_5':comp['proposed']['precision_at_k'],'selection_regret':comp['proposed']['selection_regret'],
             'different_window_rows':sum(rows[i]['isoform_class']!='identical_local_windows' for i in ix),
             'context_status':'inferred_contiguous_native_insert'})
        for i in ix:
            r = rows[i]
            item = {'row_id':r['row_id'],'gene':gene,'target_pool_n':len(ix),'global_group':r['global_group'],
               'measured_activity_unclipped':r['efficacy'],'prediction':float(pred[i]),'reference_prediction':float(ref[i]),
               'rank':ranks['actual_A_plus_B'][int(i)],'reference_rank':ranks['reference'][int(i)],
               'selected_top5':ranks['actual_A_plus_B'][int(i)]<=5 if len(ix)>=5 else None,
               'context_status':'inferred_contiguous_native_insert','isoform_class':r['isoform_class']}
            for name,v in variants.items():
                item[name+'_prediction'] = float(v[i])
                item[name+'_rank'] = ranks[name][int(i)]
            item['uncertainty_contract'] = 'prior/isoform sensitivity and independent Monte Carlo draws; not a biological confidence interval'
            candidates.append(item)
    csv_out('figure3_target_candidates',candidates)
    csv_out('target_summary',per_target)
    ab = rna['comparisons']['actual_A_plus_B']
    ni = ab['group_paired_ci95'][0] >= -cfg['metrics']['noninferiority_margin_spearman']
    rank_loss = macro(ab,'ndcg_at_k','reference')-macro(ab,'ndcg_at_k')
    activity_loss = macro(ab,'selected_activity','reference')-macro(ab,'selected_activity')
    rank_ok = rank_loss <= .02 and activity_loss <= .02
    sound = all(r['sound_vs_enumerated_float_tolerance_1e_10'] is not False and not r['negative_gap_defect'] for r in ver['rows'])
    # A positive route needs every applicable gate, including useful frozen prediction and known novelty.
    decision = 'NO-GO FOR THIS SUBMISSION ROUTE' if not cap['capability_passed'] and not ver['useful_gain_claim_admissible'] else 'REQUIRES_MANUAL_DECISION'
    resources = cycle['resources']
    cap_rows = {r['arm']:r for r in cap['rows']}
    abcap = cap_rows['actual_A_plus_B']
    lines = [f'# P4 — executed actual-model decision ({cycle["started_utc"]})','',
      f'**{decision}.** All requested comparisons completed on physical GPU 1 (NVIDIA TITAN RTX), UUID `{cycle["physical_gpu_uuid"]}`. '
      f'The run performed {cycle["optimizer_steps"]:,} optimizer steps and saved trained checkpoints.','',
      'The latest request to fully run the experiments was recorded as a separate execution authorization. The historical resource ledger remains unreconciled. '
      'The earlier gated entry point and its blocked result are preserved. Before seeing these outcomes, this execution declared that it would complete all comparisons; '
      'because the capability gate failed, subsequent RNA and verifier results are diagnostics and do not satisfy the original capability-gated positive route.','',
      '| Decision question | Executed evidence |','|---|---|',
      f'| Actual architecture and nonadditivity? | Actual `CertifiedModel`, float64; no readout rescaling. Exhaustive trained-model additive-fit residual max {cap["trained_model_nonadditivity"]["max_abs"]:.6g}. |',
      f'| Independent nonlinear capacity control? | **Fail:** training R² {abcap["train"]["r2"]:.4f}, required ≥0.95. Test MSE {abcap["test"]["mse"]:.6f}; pairwise ridge {cap_rows["pairwise_ridge"]["test"]["mse"]:.6f}. |',
      f'| RNA predictive noninferiority? | {"Pass" if ni else "Not established"}: Δ Spearman {ab["delta_spearman"]:+.4f}, paired group 95% interval [{ab["group_paired_ci95"][0]:+.4f}, {ab["group_paired_ci95"][1]:+.4f}], margin −0.02. |',
      f'| RNA candidate utility? | {"Pass" if rank_ok else "Fail"} practical criterion: macro NDCG loss {rank_loss:+.4f}; selected activity loss {activity_loss:+.4f}; each permitted loss ≤0.02. |',
      f'| Useful learned-versus-fixed verification gain? | {ver["useful_gain_claim_admissible"]}. Capability gate failed; equal-total-cost measurements below remain diagnostic. |',
      f'| Numerical certification? | `float_unverified`; exhaustive float check at 1e−10 {"found no violations" if sound else "found defects"}. This is not a numerical-enclosure proof. |',
      '| Context and chemistry? | Conditional inferred shared contiguous native insert; label-independent transcript rule. No verified construct, no chemistry-transfer test. |',
      '| Defensible ML operator distinction? | Exact support/count/partition over a declared noncrossing family plus a nonlinear residual. Edge conditioning already appears in Hojny 2024 §3.4; this experiment earns no new operator or learned-bound claim. |',
      f'| Resources and remaining? | {resources["wall_s"]:.1f} s wall, {resources["cpu_core_hours"]:.5f} measured CPU core-hours including waited children, {resources["gpu_hours"]:.5f} conservatively charged GPU-hours. Historical remaining balance unknown. |','',
      '## Independent capability result','',
      'One prespecified nonlinear motif task, 96 ensemble observations / 48 independent paired groups, exhaustive valid-member expectations. '
      'Equal-marginal pairs differ in higher-order statistics and labels. The loss is applied after averaging structure predictions. '
      'One fixed seed; 200 anchor epochs, 400 residual epochs; no outcome-driven rerun.','',
      '| Arm | Train R² | Held-out MSE |','|---|---:|---:|']
    for arm,r in cap_rows.items(): lines.append(f'| {arm} | {fmt(r["train"]["r2"])} | {fmt(r["test"]["mse"],6)} |')
    lines += ['',f'Development selected γ={cap["chosen_gamma"]}. The fit failure is concrete evidence about this bounded optimization, not proof that the architecture can never learn the task.','',
      '## RNA prediction and ranking','',
      f'1,130 matched conditional reporter rows; train/development/test = {rna["split"]["sizes"]}; global shared-13mer component split, leakage 0. '
      f'Exploratory evaluation: this dataset was already inspected in P3. Reference selected by development MSE: `{rna["reference_selected_on_development"]}`. '
      f'Actual A and B each trained for 20 epochs; development selected γ={rna["gamma"]}. All labels remain unclipped; assay identity is YFP / H1299 / 48 h.','',
      '| Arm | Spearman | Δ vs reference | Macro target NDCG@5 | Selected activity |','|---|---:|---:|---:|---:|']
    for r in predictive: lines.append(f'| {r["arm"]} | {fmt(r["spearman"])} | {r["delta_spearman"]:+.4f} | {fmt(r["macro_target_ndcg_at_5"])} | {fmt(r["macro_selected_activity"])} |')
    lines += ['', 'Ranking averages include only target pools with at least five held-out candidates. Target counts, prevalence, P@5, regret, '
              'paired group uncertainty, target-specific comparisons and leave-one-target-out sensitivity are retained in the CSV/JSON evidence. '
              'This does not establish performance on unseen targets.','',
              '## Prior, isoform and sampling sensitivity','',
              '| Declared temperature | Max prediction change | Δ held-out Spearman vs T=1 | Mean prediction Monte Carlo SE |',
              '|---|---:|---:|---:|']
    for t,r in rna['prior_sensitivity'].items(): lines.append(f'| {t} | {r["max_prediction_change"]:.5f} | {fmt(r["ranking"]["delta_spearman"])} | {r["mean_prediction_se"]:.6f} |')
    ss = rna['sampling']
    lines += ['',f'The declared prior is θ=0.5/T, without thermodynamic or expression interpretation. For the frozen γ=1 residual, mean sampling SE '
      f'is {ss["mean_raw_gamma1_se_16"]:.6f} at 16 draws and {ss["mean_raw_gamma1_se_64"]:.6f} at 64. '
      f'At the selected γ, the independent 64-draw replicate changes a prediction by at most {ss["independent_64_draw_max_prediction_change"]:.6f}. '
      f'Alternative local windows change held-out Spearman by {fmt(rna["isoform_alternative"]["delta_spearman"])}. '
      'These SEs measure inhibition-fraction prediction sampling error, not biological uncertainty. Training uses eight independently drawn valid graphs per observation and averages before loss; '
      'gradient sampling variance was not estimated in the training run. The candidate dataset records within-target ranks under each prior, alternative window and evaluation draw.','',
      '## Frozen-predictor verifier','',
      'Every method uses the same saved capability checkpoint. The original and wider sequence strata were fixed in the configuration. '
      'Already-decided cases use the skip action. Setup, policy inference, oracle calls and completed updates are timed per arm. '
      'Policy teacher generation and fitting are charged separately and amortized; primary comparison is at 100 instances. '
      'Tiny-family exhaustive evaluation is the strong applicable external comparator. CPU bounds use the frozen GPU-trained weights; no second GPU is used.','',
      f'Offline policy work: {ver["offline_teacher_and_fit_wall_s"]:.6f} s. Coverage on the original family, including its allocated offline charge:','',
      '| Total seconds | Combined | Fixed conditioning | Deterministic | Learned | Exhaustive |','|---:|---:|---:|---:|---:|---:|']
    for r in ver['equal_total_cost']:
        if r['amortization_instances']==100:
            lines.append('| '+f'{r["budget_s"]:.2f}'+' | '+' | '.join(f'{r["coverage"][a]:.3f}' for a in cfg['verifier']['arms'])+' |')
    lines += ['', 'The JSON retains both strata, exact float optima, unavailable gaps, tail times and timeouts with denominators. '
              'Neither ensemble-mean prediction nor a model threshold bound certifies measured biological ordering. Unsupported joint RNA uncertainty uses conservative independent families.','',
      '## Provenance, reproduction and remaining blocker','',
      '[The preceding audit and correction note](p4_actual_model_decision-20260911-064446-098490.md) contains the P1–P3 capability table, '
      'claims-to-evidence corrections, row exclusions and legitimate Davis acquisition attempts. Historical records and actual model code were preserved; '
      'all source hashes captured before this run were unchanged on completion. GPU batching was checked against sequential actual-model outputs and gradients; '
      'exact prior marginals and seeded samples were checked against the original oracle/sampler. Checkpoints and each stage have SHA256 sidecars.','',
      'Implemented commands, from the repository root:','',
      '```bash','python -u -m scmp.revision.p4_gpu_run',f'python -m scmp.revision.p4_gpu_report {cycle_path}','```','',
      f'Run record: [`{Path(cycle_path).name}`](../{cycle_path}). Config: [`revision_p4_gpu_v1.yaml`](../configs/revision_p4_gpu_v1.yaml).','']
    for name,item in exported.items(): lines.append(f'- [{name}](../{item["path"]}) — {item["rows"]} rows, SHA256 sidecar.')
    lines += ['',f'Peak CUDA allocation: {resources["cuda_peak_allocated_bytes"]/2**20:.1f} MiB; memory fraction cap 15%; one CPU affinity, BLAS/OpenMP thread count 1, no training workers. '
      'The execution had an automatic 900-second CPU/wall cap. Remaining local cap is not a reconciled historical allowance; no further fitting was launched.','',
      'Davis remains unavailable: `gkaf479_supplemental_files.zip` (Supplementary Table S1), Oxford Academic NAR article → Supplementary data, '
      'destination `data/davis2025/gkaf479_supplemental_files.zip`. No challenge was bypassed. Chemistry transfer, native-versus-reporter comparisons and missing dose/time/replicate covariates remain blocked.','',
      '**Recommended next action:** pivot this submission route; preserve these results as the stopping evidence rather than launch another search for favorable seeds or hyperparameters.','']
    report = Path('docs') / ('p4_actual_model_executed-'+cycle['started_utc']+'.md')
    with report.open('x') as f: f.write('\n'.join(lines))
    Path(str(report)+'.sha256').write_text(f'{sha(report)}  {report.name}\n')
    manifest = {'cycle':cycle_path,'decision':decision,'report':str(report),'report_sha256':sha(report),'exports':exported,
                'rna_spearman_noninferiority':ni,'rna_ranking_criterion':rank_ok,'numerical_float_audit_no_defects':sound,
                'resources':usage_delta(start),'training_executed':False}
    path = str(outprefix)+'-report.json'
    write_record(path,manifest)
    print(str(report)); print(path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('cycle')
    main(parser.parse_args().cycle)
