"""P4: immutable audit, measured no-fit dry run, resource-gated cycle."""
from __future__ import annotations
import argparse
import datetime
import json
import os
from pathlib import Path

# Set before numerical imports. No GPU, subprocesses or data-loader workers.
for _key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
             'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS'):
    os.environ[_key] = '1'
os.environ['CUDA_VISIBLE_DEVICES'] = ''
if hasattr(os, 'sched_getaffinity'):
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})

from .p4_resources import Budget, BudgetStop, inventory, sha, usage, usage_delta, write_record


def main(argv=None):
    ap = argparse.ArgumentParser(prog='python -m scmp.revision.p4_decision')
    ap.add_argument('--config', default='configs/revision_p4.yaml')
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument('--audit', action='store_true', help='read-only reconciliation; no model fits')
    mode.add_argument('--dry-run', action='store_true', help='measured actual-model forward/backward; zero optimizer steps')
    mode.add_argument('--run', action='store_true', help='requires verified remaining resources; refuses before fitting otherwise')
    ap.add_argument('--out-dir', default='runs/revision')
    a = ap.parse_args(argv)
    import yaml
    start = usage()
    cfg = yaml.safe_load(Path(a.config).read_text())
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%d-%H%M%S-%f')
    tag = 'audit' if a.audit else 'dryrun' if a.dry_run else 'cycle'
    prefix = Path(a.out_dir)/f'p4_{tag}-{stamp}'
    protected = [p for root in ('scmp', 'certmp', 'experiments', 'data/models')
                 for p in Path(root).rglob('*') if p.is_file() and '__pycache__' not in str(p)
                 and p.suffix in ('.py', '.json', '.pt', '.pth') and not p.name.startswith('p4_')]
    before = {str(p): sha(p) for p in protected}
    ledger = inventory(cfg)
    res = {'kind': f'p4_{tag}', 'created_utc': stamp, 'config': cfg, 'config_path': a.config,
           'config_sha256': sha(a.config), 'resource_ledger': ledger,
           'source_hashes': {str(p): sha(p) for p in Path('scmp/revision').glob('p4_*.py')},
           'data_hashes': {str(p): sha(p) for p in Path('data').rglob('*') if p.is_file()},
           'protected_hashes_before': before, 'training_executed': False,
           'numerical_status': 'float_unverified'}
    code = 0
    if a.run and not ledger['admitted']:
        res.update(status='BLOCKED_RESOURCE_RECONCILIATION',
                   blocked_run='independent ensemble capability training; downstream RNA and learned verifier gated',
                   decision='NOT_EVALUATED',
                   decision_reason='an unexecuted comparison cannot be called GO or failed')
        code = 2
    else:
        import torch
        torch.set_num_threads(1)
        torch.set_num_interop_threads(1)
        limit = cfg['resource']['audit_wall_seconds']
        if a.run:
            limit = min(cfg['resource']['max_cycle_cpu_seconds'], cfg['resource']['max_cycle_wall_seconds'],
                        ledger['verified_spendable_cpu_seconds'], ledger['verified_spendable_wall_seconds'])
        budget = Budget(limit, limit)
        try:
            with budget.enforce():
                if a.audit:
                    from .p4_data import historical_evidence, audit_rows
                    res['historical'] = historical_evidence(cfg)
                    data = audit_rows(cfg)
                    data_path = str(prefix)+'-rows.json'
                    write_record(data_path, data)
                    res.update(row_record=data_path, row_record_sha256=sha(data_path),
                               data_summary=data['summary'], status='AUDIT_COMPLETE_NO_FITS')
                elif a.dry_run:
                    from .p4_capability import smoke
                    res['smoke'] = smoke(cfg['capability'], budget)
                    res.update(status='DRY_RUN_COMPLETE_NO_FITS', training_admitted=ledger['admitted'])
                else:
                    expected_instances = cfg['capability']['n_instances'] + 2*(cfg['verifier']['development_instances']+cfg['verifier']['test_instances'])
                    if expected_instances > ledger['verified_spendable_instances']:
                        raise BudgetStop('insufficient reconciled instance allowance')
                    from .p4_capability import smoke
                    res['measured_preflight'] = smoke(cfg['capability'], budget)
                    if res['measured_preflight']['estimated_control_training_wall_s'] > limit/2:
                        raise BudgetStop('measured control estimate exceeds half the cycle reservation')
                    write_record(str(prefix)+'-preregistration.json', res)
                    from .p4_capability import run_control
                    res['training_executed'] = True
                    res['training_instances'] = expected_instances
                    res['control'], model = run_control(cfg['capability'], budget)
                    write_record(str(prefix)+'-control.json', res['control'])
                    if res['control']['capability_passed']:
                        from .p4_verifier import run_verifier
                        res['verifier'] = run_verifier(model, cfg, budget)
                        write_record(str(prefix)+'-verifier.json', res['verifier'])
                        from .p4_rna import run_rna
                        res['rna'] = run_rna(cfg, budget)
                        res['decision'] = 'DECISION_REQUIRES_NOVELTY_AND_UTILITY_REVIEW'
                    else:
                        res['rna'] = {'status': 'NOT_RUN_CAPABILITY_GATE_FAILED'}
                        res['verifier'] = {'status': 'NOT_RUN_NO_USEFUL_PREDICTOR'}
                        res['decision'] = 'CAPABILITY_FAILED_VERIFIER_NOT_EVALUATED'
                    res['status'] = 'BOUNDED_CYCLE_COMPLETE'
        except BudgetStop as exc:
            res.update(status='STOPPED_AT_RESOURCE_CAP', stopping_reason=str(exc), decision='NOT_EVALUATED_INCOMPLETE')
            code = 3
    after = {p: sha(p) for p in before}
    res.update(protected_hashes_after=after, preserved=before == after, resources=usage_delta(start))
    if before != after:
        raise RuntimeError('protected model/source drift detected')
    res['resources']['training_instances'] = res.get('training_instances', 0)
    path = write_record(str(prefix)+'.json', res)
    print(json.dumps({'record': path, 'status': res['status'], 'resources': res['resources'],
                      'training_executed': res['training_executed'], 'preserved': res['preserved']}, indent=2))
    return code


if __name__ == '__main__':
    raise SystemExit(main())
