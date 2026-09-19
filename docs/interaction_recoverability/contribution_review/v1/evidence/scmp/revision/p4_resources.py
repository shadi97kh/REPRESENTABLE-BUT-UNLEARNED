"""P4 accounting. Missing historical measurements never become free compute."""
from __future__ import annotations

import contextlib
import hashlib
import json
import os
from pathlib import Path
import resource
import signal
import time


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def checked_record(path):
    path = Path(path)
    side = Path(str(path) + '.sha256')
    if not side.exists() or side.read_text().split()[0] != sha(path):
        raise ValueError(f'missing or incorrect SHA256 sidecar: {path}')
    return json.loads(path.read_text())


def write_record(path, obj):
    """Exclusive creation: never overwrite a result, including within one second."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as f:
        json.dump(obj, f, indent=2, sort_keys=True, allow_nan=False)
        f.write('\n')
    with Path(str(path) + '.sha256').open('x') as f:
        f.write(f'{sha(path)}  {path.name}\n')
    return str(path)


def usage():
    a, b = resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)
    return dict(wall_s=time.monotonic(), process_cpu_s=a.ru_utime + a.ru_stime,
                child_cpu_s=b.ru_utime + b.ru_stime, max_rss_kib=a.ru_maxrss)


def usage_delta(start):
    end = usage()
    out = {k: end[k] - start[k] for k in ('wall_s', 'process_cpu_s', 'child_cpu_s')}
    out.update(max_rss_kib=end['max_rss_kib'], gpu_hours=0, child_processes_spawned=0,
               torch_threads=1, blas_threads=1, openmp_threads=1)
    out['cpu_core_hours'] = (out['process_cpu_s'] + out['child_cpu_s']) / 3600
    return out


class BudgetStop(RuntimeError):
    pass


class Budget:
    """One pinned CPU, no workers/GPU. Wall alarm also interrupts long native calls.

    CPU is measured with getrusage (including waited children); the sole CPU
    affinity and no-child execution contract put CPU use below the wall cap.
    Calls into training additionally check actual process+child CPU each step.
    """
    def __init__(self, cpu_s, wall_s):
        self.cpu_s, self.wall_s = float(cpu_s), float(wall_s)
        self.start = usage()

    def check(self):
        u = usage_delta(self.start)
        if u['process_cpu_s'] + u['child_cpu_s'] >= self.cpu_s or u['wall_s'] >= self.wall_s:
            raise BudgetStop('P4 automatic resource cap reached')

    @contextlib.contextmanager
    def enforce(self):
        def stop(*_):
            raise BudgetStop('P4 wall timer stopped the current stage')
        old_handler = signal.signal(signal.SIGALRM, stop)
        previous = signal.getitimer(signal.ITIMER_REAL)
        duration = min(self.cpu_s, self.wall_s)
        if previous[0]:
            duration = min(duration, previous[0])
        signal.setitimer(signal.ITIMER_REAL, duration)
        began = time.monotonic()
        try:
            yield self
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
            signal.signal(signal.SIGALRM, old_handler)
            if previous[0]:
                signal.setitimer(signal.ITIMER_REAL,
                                 max(0.000001, previous[0] - (time.monotonic() - began)), previous[1])


def inventory(cfg):
    import yaml
    auth_path = cfg['resource']['authorization']
    auth = yaml.safe_load(Path(auth_path).read_text())['pilot']['authorized_budget']
    paths, rows, unknown, p4_usage = [], [], [], []
    # All historical run records are inventoried, not just the latest successful run.
    for p in sorted(Path('runs').rglob('*.json')):
        if p.name.startswith('p4_') or 'p4-' in str(p):
            prior = json.loads(p.read_text())
            if isinstance(prior, dict) and 'resources' in prior:
                p4_usage.append({'path': str(p), 'sha256': sha(p), **prior['resources']})
            continue
        paths.append({'path': str(p), 'sha256': sha(p)})
        d = json.loads(p.read_text())
        vals = d if isinstance(d, list) else [d]
        times = [v for v in vals if isinstance(v, dict) and
                 any(k in v for k in ('runtime_s', 'elapsed_s', 'cpu_seconds_total'))]
        if not times:
            unknown.append(str(p))
        for v in times:
            if v.get('reparsed_from'):
                # A reparse references the same execution, not another execution.
                continue
            wall = v.get('runtime_s', v.get('elapsed_s', v.get('cpu_seconds_total')))
            rows.append({'path': str(p), 'wall_s_or_legacy_cpu_proxy': wall,
                         'reported_cpu_core_hours': v.get('cpu_core_hours', v.get('core_hours')),
                         'reported_gpu_hours': v.get('gpu_hours'),
                         'actual_process_cpu_s': v.get('process_cpu_s'),
                         'actual_child_cpu_s': v.get('child_cpu_s'),
                         'legacy_cpu_seconds_total_is_wall_proxy': 'cpu_seconds_total' in v,
                         'thread_accounting_verified': False})
    # A constructive counterexample to treating reported wall-hours as CPU-hours:
    # just P3 plus BOTH approved original pilots exceed the CPU cap at 16 threads.
    selected = [r for r in rows if r['path'].startswith('runs/pilot/pilot-') or
                r['path'].endswith('rna_pilot-20260906-011411.json')]
    upper_example = sum(r['wall_s_or_legacy_cpu_proxy'] * 16 / 3600 for r in selected)
    errors = ['historical process/child CPU and BLAS/OpenMP counts are not recorded',
              'GPU runs report zero CPU charge; this does not measure host work',
              'wall/instance cap scope cannot be extended from per-run checks',
              'some executions, failures, and tests have no complete resource record']
    result = {'authorization_path': auth_path, 'authorization_sha256': sha(auth_path),
              'authorization': auth, 'historical_files': paths, 'recorded_executions': rows,
              'p4_usage_records': p4_usage,
              'records_without_top_level_runtime': unknown,
              'reported_wall_proxy_hours': sum(r['wall_s_or_legacy_cpu_proxy'] for r in rows) / 3600,
              'reported_gpu_hours_revision_only': sum(r['reported_gpu_hours'] or 0 for r in rows),
              'counterexample_16_core_hours_pilots_plus_p3': upper_example,
              'counterexample_records': [r['path'] for r in selected],
              'errors': errors, 'balance_status': 'unverified_not_proven_exhausted',
              'verified_spendable_cpu_seconds': 0.0, 'admitted': False,
              'note': 'User confirmed no reconciled ledger on 2026-09-11. No new cap is inferred.'}
    ledger_path = cfg['resource'].get('reconciliation')
    if ledger_path:
        ledger = checked_record(ledger_path)
        if ledger.get('authorization_sha256') != sha(auth_path):
            raise ValueError('reconciliation must reference the existing authorization')
        if ledger.get('historical_files') != paths:
            raise ValueError('reconciliation is stale or omits historical files')
        if not all(ledger.get('accounting_verified', {}).get(k) is True
                   for k in cfg['resource']['required_accounting']):
            raise ValueError('incomplete resource reconciliation')
        # Each charge must be supported by a local, checksum-checked measurement.
        if not ledger.get('measurement_records'):
            raise ValueError('reconciliation has no supporting measurements')
        for rec in ledger['measurement_records']:
            checked_record(rec)
        consumed = ledger['consumed']
        remaining = {k: auth[k] - consumed[k] for k in
                     ('cpu_core_hours', 'gpu_hours', 'wall_clock_hours', 'max_instances')}
        remaining['cpu_core_hours'] -= sum(r['cpu_core_hours'] for r in p4_usage)
        remaining['wall_clock_hours'] -= sum(r['wall_s'] for r in p4_usage)/3600
        remaining['gpu_hours'] -= sum(r['gpu_hours'] for r in p4_usage)
        remaining['max_instances'] -= sum(r.get('training_instances', 0) for r in p4_usage)
        if any(v < 0 for v in consumed.values()) or any(v < 0 for v in remaining.values()):
            raise ValueError('negative usage or exhausted existing authorization')
        fraction = min(0.5, cfg['resource']['fraction_of_verified_remaining'])
        result.update(remaining=remaining, balance_status='reconciled', errors=[],
                      verified_spendable_cpu_seconds=remaining['cpu_core_hours'] * 3600 * fraction,
                      verified_spendable_wall_seconds=remaining['wall_clock_hours'] * 3600 * fraction,
                      verified_spendable_instances=int(remaining['max_instances'] * fraction),
                      admitted=remaining['cpu_core_hours'] > 0 and remaining['wall_clock_hours'] > 0)
    return result
