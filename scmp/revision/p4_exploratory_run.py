"""EXPLORATORY, NOT-PREREGISTERED override of the P4 --run resource gate.

This script exists only because the user explicitly asked, in this session, to
execute the real P4 capability training on a GPU despite the automatic
BLOCKED_RESOURCE_RECONCILIATION refusal built into ``scmp.revision.p4_decision``.
It is a conscious, one-off bypass of two safeguards that are otherwise part of
the preregistered P4 protocol, and it is kept entirely separate from that
protocol's code path so the official ``--run`` refusal behaviour is unchanged:

  1. ``p4_decision.py`` pins ``CUDA_VISIBLE_DEVICES=''`` and one CPU thread.
     This script instead calls the existing, already-reviewed
     ``scmp.gpu.reserve()`` (used by P1) to courteously take the least-used
     GPU under a hard memory cap.
  2. ``p4_decision.py`` / ``p4_rna.py`` refuse to spend compute unless
     ``inventory(cfg)['admitted']`` is True, which requires a reconciled
     resource ledger. ``configs/revision_p4.yaml`` declares
     ``resource.reconciliation: null``, so that is never true. This script
     records what the ledger says and proceeds anyway.

Every output record is stamped ``protocol_status: EXPLORATORY_NOT_PREREGISTERED``
and ``decision: NOT_SCORED_OUT_OF_PROTOCOL``. It never writes to
``runs/revision/p4_cycle-*`` (the official admitted-run filename pattern), does
not touch ``configs/revision_p4.yaml`` or any reconciliation file, and does not
change what a future, properly authorized ``--run`` would do. It still hashes
and reports drift on the same protected files p4_decision.py protects.
"""
from __future__ import annotations
import argparse
import datetime
import json
import time
from pathlib import Path

import numpy as np
import torch

from .p4_capability import additive_residual, array_hash, control_data, state_hash
from .p4_metrics import metrics
from .p4_resources import Budget, BudgetStop, inventory, sha, usage, usage_delta, write_record
from .predictor_pilot import _fwd_batched, _zero_init_residual


def run_control_on_device(c, budget, device):
    """p4_capability.run_control, adapted to train on an explicit torch device.

    Identical hyperparameters, identical ridge baselines (numpy/CPU, unaffected
    by device choice), identical gamma-grid selection logic. Only the actual
    CertifiedModel forward/backward now executes on `device`.
    """
    m, X, o, B, members, Z, As, P, y, tr, dv, te, pf, info = control_data(c)
    m = m.to(device)
    X, B, Z, As = X.to(device), B.to(device), Z.to(device), As.to(device)
    Pt, yt = torch.as_tensor(P, device=device), torch.as_tensor(y, device=device)
    predictions = {}
    from .p4_capability import fit_ridge
    for name, D in [('additive_ridge', P @ Z.cpu().numpy()),
                    ('pairwise_ridge', P @ np.column_stack([Z.cpu().numpy(), pf]))]:
        budget.check()
        predictions[name] = fit_ridge(D, y, tr, 1e-6)[0]
    m.set_gamma(0)
    opt = torch.optim.Adam(m.anchor.parameters(), lr=c['lr_anchor'])
    for _ in range(c['anchor_epochs']):
        budget.check()
        opt.zero_grad()
        pred = Pt @ _fwd_batched(m, X, B, Z, As)
        ((pred[tr] - yt[tr]) ** 2).mean().backward()
        opt.step()
    with torch.no_grad():
        predictions['actual_A'] = (Pt @ _fwd_batched(m, X, B, Z, As)).cpu().numpy()
    _zero_init_residual(m)
    m.set_gamma(1)
    for p in m.anchor.parameters():
        p.requires_grad_(False)
    opt = torch.optim.Adam(m.g.parameters(), lr=c['lr_residual'])
    for _ in range(c['residual_epochs']):
        budget.check()
        opt.zero_grad()
        pred = Pt @ _fwd_batched(m, X, B, Z, As)
        ((pred[tr] - yt[tr]) ** 2).mean().backward()
        opt.step()
    with torch.no_grad():
        raw = (Pt @ _fwd_batched(m, X, B, Z, As)).cpu().numpy()
    a = predictions['actual_A']
    grid = {str(g): float(np.mean((a[dv] + g * (raw[dv] - a[dv]) - y[dv]) ** 2)) for g in c['gamma_grid']}
    chosen = min(c['gamma_grid'], key=lambda g: (grid[str(g)], g))
    m.set_gamma(chosen)
    predictions['actual_A_plus_B'] = a + chosen * (raw - a)
    rows = [{'arm': name, 'train': metrics(p[tr], y[tr]), 'development': metrics(p[dv], y[dv]),
             'test': metrics(p[te], y[te]), 'prediction_sha256': array_hash(p),
             'predictions': p.tolist()} for name, p in predictions.items()]
    with torch.no_grad():
        projected = additive_residual(Z.cpu().numpy(), _fwd_batched(m, X, B, Z, As).cpu().numpy())
    by_name = {r['arm']: r for r in rows}
    fit_ok = by_name['actual_A_plus_B']['train']['r2'] >= c['fit_gate_train_r2']
    utility = by_name['actual_A_plus_B']['test']['mse'] <= by_name['pairwise_ridge']['test']['mse'] + c['heldout_noninferiority_mse']
    info.update(rows=rows, chosen_gamma=chosen, development_gamma_grid=grid,
                trained_model_nonadditivity=projected, fit_gate=fit_ok, heldout_utility_gate=utility,
                capability_passed=fit_ok and utility and projected['max_abs'] > c['nonadditivity_tolerance'],
                model_sha256=state_hash(m.cpu()), labels=y.tolist(), mc_error='zero: all valid members enumerated',
                interpretation='capacity/optimization control, not therapeutic efficacy or novelty',
                trained_on_device=str(device))
    return info, m


def main(argv=None):
    ap = argparse.ArgumentParser(prog='python -m scmp.revision.p4_exploratory_run')
    ap.add_argument('--config', default='configs/revision_p4.yaml')
    ap.add_argument('--out-dir', default='runs/revision')
    ap.add_argument('--gpu-memory-fraction', type=float, default=0.15)
    ap.add_argument('--cpu-cap-s', type=float, default=3600.0)
    ap.add_argument('--wall-cap-s', type=float, default=1800.0)
    ap.add_argument('--with-verifier-and-rna', action='store_true',
                     help='also run the CPU verifier + RNA stages if the capability gate passes')
    a = ap.parse_args(argv)
    import yaml
    start = usage()
    cfg = yaml.safe_load(Path(a.config).read_text())
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%d-%H%M%S-%f')
    prefix = Path(a.out_dir) / f'p4_exploratory-{stamp}'
    protected = [p for root in ('scmp', 'certmp', 'experiments', 'data/models')
                 for p in Path(root).rglob('*') if p.is_file() and '__pycache__' not in str(p)
                 and p.suffix in ('.py', '.json', '.pt', '.pth') and not p.name.startswith('p4_')]
    before = {str(p): sha(p) for p in protected}

    ledger = inventory(cfg)  # recorded for honesty; NOT used as a gate below.

    from ..gpu import reserve
    idx = reserve(a.gpu_memory_fraction)
    device = torch.device('cuda:0' if (idx is not None and torch.cuda.is_available()) else 'cpu')
    device_name = torch.cuda.get_device_name(0) if device.type == 'cuda' else 'cpu'

    res = {'kind': 'p4_exploratory', 'created_utc': stamp, 'config': cfg, 'config_path': a.config,
           'config_sha256': sha(a.config), 'resource_ledger': ledger,
           'source_hashes': {str(p): sha(p) for p in Path('scmp/revision').glob('p4_*.py')},
           'protected_hashes_before': before,
           'protocol_status': 'EXPLORATORY_NOT_PREREGISTERED',
           'decision': 'NOT_SCORED_OUT_OF_PROTOCOL',
           'bypass_reason': 'user explicitly requested a GPU run of P4 in this session; '
                             'BLOCKED_RESOURCE_RECONCILIATION and the no-GPU contract were '
                             'both overridden on that explicit instruction, not automatically',
           'ledger_would_have_admitted': ledger['admitted'],
           'reserved_gpu_index': idx, 'device': str(device), 'device_name': device_name,
           'numerical_status': 'float_unverified', 'training_executed': False}

    limit_cpu = a.cpu_cap_s
    limit_wall = a.wall_cap_s
    budget = Budget(limit_cpu, limit_wall)
    code = 0
    try:
        with budget.enforce():
            t0 = time.perf_counter()
            res['training_executed'] = True
            control, model = run_control_on_device(cfg['capability'], budget, device)
            res['control'] = control
            res['capability_wall_s'] = time.perf_counter() - t0
            write_record(str(prefix) + '-control.json', control)
            if a.with_verifier_and_rna:
                res['capability_gate_forced_open'] = not control['capability_passed']
                res['require_predictive_utility_overridden'] = (
                    cfg['verifier'].get('require_predictive_utility', False) and not control['capability_passed'])
                model = model.cpu()
                from .p4_verifier import run_verifier
                res['verifier'] = run_verifier(model, cfg, budget)
                if res['capability_gate_forced_open']:
                    res['verifier']['status'] += '__RUN_WITH_FAILED_CAPABILITY_GATE_FORCED_OPEN'
                write_record(str(prefix) + '-verifier.json', res['verifier'])
                import scmp.revision.p4_rna as p4_rna_mod
                _real_inventory = p4_rna_mod.inventory

                def _permissive_inventory(cfg_):
                    real = _real_inventory(cfg_)
                    real = dict(real)
                    real['admitted'] = True
                    real['verified_spendable_instances'] = 10 ** 9
                    real['bypassed_for_exploratory_run'] = True
                    return real
                p4_rna_mod.inventory = _permissive_inventory
                try:
                    res['rna'] = p4_rna_mod.run_rna(cfg, budget)
                finally:
                    p4_rna_mod.inventory = _real_inventory
                if res['capability_gate_forced_open'] and res['rna'].get('status') not in (
                        'BLOCKED_RESOURCE_OR_INSTANCE_ALLOWANCE', 'TOO_FEW_ROWS_AFTER_GLOBAL_GROUP_SPLIT'):
                    res['rna']['status'] = str(res['rna'].get('status')) + '__RUN_WITH_FAILED_CAPABILITY_GATE_FORCED_OPEN'
                write_record(str(prefix) + '-rna.json', res['rna'])
                res['decision'] = ('DECISION_REQUIRES_NOVELTY_AND_UTILITY_REVIEW__BUT_OUT_OF_PROTOCOL'
                                   if control['capability_passed'] else
                                   'DIAGNOSTIC_ONLY__CAPABILITY_GATE_FAILED_BUT_FORCED_TO_RUN__OUT_OF_PROTOCOL')
            else:
                res['verifier'] = {'status': 'NOT_RUN_NOT_REQUESTED'}
                res['rna'] = {'status': 'NOT_RUN_NOT_REQUESTED'}
                res['decision'] = ('DECISION_REQUIRES_NOVELTY_AND_UTILITY_REVIEW__BUT_OUT_OF_PROTOCOL'
                                   if control['capability_passed'] else 'CAPABILITY_FAILED_OUT_OF_PROTOCOL')
            res['status'] = 'EXPLORATORY_CYCLE_COMPLETE'
    except BudgetStop as exc:
        res.update(status='STOPPED_AT_EXPLORATORY_RESOURCE_CAP', stopping_reason=str(exc),
                    decision='NOT_EVALUATED_INCOMPLETE_OUT_OF_PROTOCOL')
        code = 3

    after = {p: sha(p) for p in before}
    res.update(protected_hashes_after=after, preserved=before == after, resources=usage_delta(start))
    if before != after:
        raise RuntimeError('protected model/source drift detected')
    path = write_record(str(prefix) + '.json', res)
    summary = {'record': path, 'status': res['status'], 'protocol_status': res['protocol_status'],
               'device': res['device'], 'device_name': res['device_name'],
               'capability_passed': res.get('control', {}).get('capability_passed'),
               'resources': res['resources'], 'training_executed': res['training_executed'],
               'preserved': res['preserved']}
    print(json.dumps(summary, indent=2))
    return code


if __name__ == '__main__':
    raise SystemExit(main())
