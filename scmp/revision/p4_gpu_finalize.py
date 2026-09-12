"""Finalize immutable executed evidence after no-fit checks; no experimental rerun."""
import argparse
import csv
from pathlib import Path
import numpy as np
from .p4_resources import checked_record,sha,write_record,usage,usage_delta
from .p4_verifier import decision


def main(cycle_path, test_path):
    start=usage()
    cycle=checked_record(cycle_path)
    prefix=cycle_path.replace('-cycle.json','')
    post=checked_record(prefix+'-postcheck.json')
    exported=checked_record(prefix+'-report.json')
    tests=checked_record(test_path)
    assert tests['exit_code']==0 and post['frozen_checkpoint_unchanged']
    ver=checked_record(cycle['artifacts']['verifier']['path'])
    cfg=checked_record(cycle['artifacts']['prerun']['path'])['config']
    # Recompute coverage AND gaps from the same accepted event under total cost.
    data=[]
    for stratum in ['original_family','wider_family']:
        for cap in cfg['verifier']['budgets_s']:
            for n in cfg['verifier']['offline_amortization_instances']:
                for arm in cfg['verifier']['arms']:
                    rr=[r for r in ver['rows'] if r['case']['stratum']==stratum and r['cap_s']==cap and r['arm']==arm]
                    charge=ver['offline_teacher_and_fit_wall_s']/n if arm=='learned_selective' else 0.
                    resolved=0; gaps=[]; unresolved_gaps=[]; missing=0
                    for r in rr:
                        events=[e for e in r['events'] if e['time_s']<=cap-charge]
                        if not events:
                            missing+=1; continue
                        e=events[-1]; gap=e['hi']-e['lo']
                        is_decided=decision(e['lo'],e['hi'],r['case']['tau'])
                        resolved+=is_decided; gaps.append(gap)
                        if not is_decided: unresolved_gaps.append(gap)
                    data.append({'stratum':stratum,'arm':arm,'total_budget_s':cap,'amortization_instances':n,
                      'offline_charge_s':charge,'n':len(rr),'decision_coverage':resolved/len(rr),
                      'mean_gap_at_total_budget':float(np.mean(gaps)) if gaps else None,
                      'gap_unavailable_n':missing,'unresolved_n':len(rr)-resolved,
                      'mean_unresolved_gap_when_available':float(np.mean(unresolved_gaps)) if unresolved_gaps else None,
                      'unresolved_gap_available_n':len(unresolved_gaps),
                      'online_timeouts':sum(r['timeout'] for r in rr),
                      'online_wall_p90_s':float(np.quantile([r['wall_s'] for r in rr],.9)),
                      'numerical_status':'float_unverified','capability_gate_passed':cycle['capability_passed']})
    cp=Path(prefix+'-figure2_coverage_total_compute_v2.csv')
    with cp.open('x',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(data[0]));writer.writeheader();writer.writerows(data)
    Path(str(cp)+'.sha256').write_text(f'{sha(cp)}  {cp.name}\n')
    cand_path=exported['exports']['figure3_target_candidates']['path']
    with open(cand_path) as f: candidates=list(csv.DictReader(f))
    genes=sorted({r['gene'] for r in candidates if int(r['target_pool_n'])>=5})
    changes={}
    for variant in ['prior_T_0.5','prior_T_2.0','alternative_isoform','independent_64_draw','draw_16']:
        changes[variant]=0
        for gene in genes:
            rr=[r for r in candidates if r['gene']==gene]
            a={r['row_id'] for r in rr if int(r['rank'])<=5}
            b={r['row_id'] for r in rr if int(r[variant+'_rank'])<=5}
            changes[variant]+=a!=b
    old=Path(exported['report'])
    assert sha(old)==exported['report_sha256']
    text=old.read_text().replace('gradient sampling variance was not estimated in the training run.',
        'a separate frozen-checkpoint diagnostic below quantifies conditional gradient sampling variability without further fitting.')
    text=text.replace(exported['exports']['figure2_coverage_total_compute']['path'],str(cp)).replace('— 45 rows, SHA256 sidecar.','— 90 rows, both strata and gaps at equal total cost, SHA256 sidecar.')
    text=text.replace('Implemented commands, from the repository root:',
        'Commands used, from the repository root (outputs are created exclusively; rerunning an exporter with the same output prefix refuses overwrites):')
    appendix=['','## Final checkpoint and evidence checks','',
      f'All 54 model/P4 contract tests passed ([record](../{test_path})). Three trained checkpoints reload with identical parameter hashes and finite parameters. '
      f'On actual 50-nt windows, batched outputs match sequential `CertifiedModel` outputs to {post["batched_50nt_output_max_error"]:.2g}; '
      f'exact cached marginals match the original implementation to {post["exact_50nt_marginal_max_error"]:.2g}; '
      f'reloaded checkpoints reproduce saved held-out predictions to {post["saved_prediction_reproduction_max_error"]:.2g}.','',
      'The post-run gradient diagnostic fixes the final weights, uses the first 16 training rows and four independent sample replicates, '
      'and measures the residual at its training scale γ=1. No optimizer runs and the saved state is unchanged. '
      f'Relative gradient RMS standard deviation is {post["gradient_sampling"]["8"]["relative_rms_sd"]:.2%} at 8 draws and '
      f'{post["gradient_sampling"]["32"]["relative_rms_sd"]:.2%} at 32 draws. This measures conditional sampling variability, not variation across optimization runs or biological experiments. '
      f'[Postcheck record](../{prefix}-postcheck.json).','',
      f'Top-five shortlist membership changes across {len(genes)} eligible held-out target pools:','',
      '| Perturbation | Target shortlists changed |','|---|---:|']
    appendix += [f'| {k} | {v}/{len(genes)} |' for k,v in changes.items()]
    appendix += ['','The original verifier stratum is saturated once setup completes: every method has coverage 1 at 0.05 and 0.20 seconds. '
      'This is a saturation finding, not evidence for a learned bounding layer. The prespecified wider stratum is retained below and in the expanded figure dataset; '
      'its predictive usefulness was never established. Both strata therefore remain diagnostic.','',
      '| Wider stratum total seconds | Combined | Fixed conditioning | Deterministic | Learned | Exhaustive |',
      '|---:|---:|---:|---:|---:|---:|']
    for cap in cfg['verifier']['budgets_s']:
        selected=[r for r in data if r['stratum']=='wider_family' and r['total_budget_s']==cap and r['amortization_instances']==100]
        appendix.append('| '+str(cap)+' | '+' | '.join(f'{next(r["decision_coverage"] for r in selected if r["arm"]==a):.3f}' for a in cfg['verifier']['arms'])+' |')
    records=[cycle,post,tests,exported]
    first_validation=checked_record('runs/revision/p4_gpu_validation-1789146455492715666.json')
    records.append(first_validation)
    totals={k:sum(r['resources'][k] for r in records) for k in ['process_cpu_s','child_cpu_s','cpu_core_hours','gpu_hours','wall_s']}
    appendix += ['',f'Instrumented execution and validation totals before this final export: {totals["wall_s"]:.2f} seconds, '
      f'{totals["cpu_core_hours"]:.5f} CPU core-hours and {totals["gpu_hours"]:.5f} conservatively charged GPU-hours. '
      'These measurements cover the recorded computation regions; interpreter/module imports and uninstrumented editor/read-only shell work are excluded. '
      'They do not reconcile the historical allowance. The experiment itself used 269.28 of its 900-second local automatic wall limit; no further training is authorized by this arithmetic.','',
      'Additional implemented commands used for these checks:','','```bash',
      'python -m scmp.revision.p4_checks',f'python -m scmp.revision.p4_gpu_postcheck {cycle_path}',
      f'python -m scmp.revision.p4_gpu_finalize {cycle_path} {test_path}','```','']
    report=old.with_name(old.stem+'-final.md')
    with report.open('x') as f: f.write(text+'\n'.join(appendix))
    Path(str(report)+'.sha256').write_text(f'{sha(report)}  {report.name}\n')
    result={'report':str(report),'report_sha256':sha(report),'cycle':cycle_path,
            'expanded_figure2':{'path':str(cp),'sha256':sha(cp)},'shortlist_changes':changes,'eligible_target_pools':len(genes),
            'instrumented_execution_and_checks_totals':totals,'resources':usage_delta(start),'training_executed':False,
            'optimizer_steps':0,'source_sha256':sha(__file__)}
    write_record(prefix+'-final.json',result)
    print(str(report)); print(totals)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('cycle');parser.add_argument('tests')
    a=parser.parse_args();main(a.cycle,a.tests)
