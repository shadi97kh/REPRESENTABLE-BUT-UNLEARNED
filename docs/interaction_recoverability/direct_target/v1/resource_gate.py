"""Task cap layered onto the unchanged learned-target ledger and watchdog."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[4]
RUN = ROOT / 'runs/interaction_recoverability/learned-target-v1-20260912-164000'
HERE = Path(__file__).resolve().parent
TASK = 'direct-target-v1'
TASK_CAP = 600.0
OVERHEAD = 150.0
EXPECTED_LEDGER = '29e33d348dbf42e18d3cf4dbe3ae34496af584f523ef5e868b0b87578a8dcc4a'
sys.path.insert(0, str(ROOT))
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
sys.dont_write_bytecode = True
from interaction_recoverability.learned_target import budget


def read_verified():
    rows = budget.events(RUN)
    starts = {r['job_id'] for r in rows if r['event'] == 'start'}
    ends = {r['job_id'] for r in rows if r['event'] == 'end'}
    if starts != ends:
        raise RuntimeError('Unclosed global job')
    first = rows[0]
    if first['cpu_s_cap'] != 7200 or first['aggregate_job_wall_s_cap'] != 7200 or first['reserve_s'] != 600:
        raise RuntimeError('Unexpected global allowance')
    config = json.loads((RUN / 'frozen_config.json').read_text())
    for rel, expected in config['learner_code_sha256'].items():
        if hashlib.sha256((ROOT / rel).read_bytes()).hexdigest() != expected:
            raise RuntimeError('Frozen learner changed: ' + rel)
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage', choices=['register', 'run', 'verify'])
    p.add_argument('--seconds', type=float, default=30.)
    p.add_argument('command', nargs=argparse.REMAINDER)
    a = p.parse_args()
    lock = (HERE / 'task.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    rows = read_verified()
    markers = [r for r in rows if r.get('event') == 'continuation' and r.get('task') == TASK]
    if a.stage == 'register':
        if markers:
            raise RuntimeError('Task already registered; no reset')
        if hashlib.sha256((RUN / 'ledger.jsonl').read_bytes()).hexdigest() != EXPECTED_LEDGER:
            raise RuntimeError('Pre-continuation ledger mismatch; reconcile statically')
        used = sum(r.get('conservative_charge_s', 0.) for r in rows)
        if budget.CAP - used - budget.RESERVE < TASK_CAP:
            raise RuntimeError('Insufficient verified balance preserving reserve')
        paths = []
        for base in ['docs/interaction_recoverability', 'interaction_recoverability', 'configs', 'tests', str(RUN.relative_to(ROOT))]:
            paths.extend(p for p in (ROOT / base).rglob('*') if p.is_file() and '__pycache__' not in p.parts)
        hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths
                  if HERE not in p.parents and p.name not in ['ledger.jsonl', 'lock']}
        hashes['runs/interaction_recoverability/prep-20260912-082154/ledger.jsonl'] = hashlib.sha256((ROOT / 'runs/interaction_recoverability/prep-20260912-082154/ledger.jsonl').read_bytes()).hexdigest()
        with (HERE / 'preservation.json').open('x') as f:
            json.dump(hashes, f, indent=2)
        budget.append(RUN / 'ledger.jsonl', dict(event='continuation', task=TASK,
            additional_cpu_s_cap=TASK_CAP, additional_job_wall_s_cap=TASK_CAP,
            prior_charge_s=used, prior_through_job=max(ends['job_id'] for ends in rows if ends['event']=='end'),
            prior_ledger_sha256=EXPECTED_LEDGER, conservative_charge_s=OVERHEAD,
            measured_cpu_s=None, gpu_s_cap=0, workers=1, numerical_threads=1,
            note='Explicit direct-target continuation; no allowance reset. 150s conservative debit for discovery, static source/file reads, authoring, agent static-read overhead and final reporting outside wrapped jobs. Historical ledger read only.'))
        rows = read_verified()
        markers = [r for r in rows if r.get('event') == 'continuation' and r.get('task') == TASK]
    if len(markers) != 1:
        raise RuntimeError('Register this continuation before execution')
    marker = markers[0]
    charge = sum(r.get('conservative_charge_s', 0.) for r in rows)
    task_used = charge - marker['prior_charge_s']
    remaining = min(TASK_CAP-task_used, budget.CAP-charge-budget.RESERVE)
    if a.stage == 'run':
        if not 0 < a.seconds <= 120:
            raise ValueError('Requested job deadline must be in (0,120]')
        timeout = min(a.seconds, remaining-5.)
        if timeout <= 0:
            raise RuntimeError('Task/reserve cap exhausted')
        command = a.command[1:] if a.command[:1] == ['--'] else a.command
        if not command:
            raise ValueError('Explicit command required')
        sys.argv = ['unchanged-budget', '--run-root', str(RUN), '--phase', 'compute', '--job-seconds', str(timeout), '--', *command]
        budget.main()
    else:
        completed = [r for r in rows if r['event']=='end' and r['job_id'] > marker['prior_through_job']]
        summary = dict(task=TASK, through_job=max(r['job_id'] for r in rows if r['event']=='end'),
            task_measured_cpu_s=sum(r['measured_cpu_s'] for r in completed),
            task_job_wall_s=sum(r['wall_s'] for r in completed), task_conservative_charge_s=task_used,
            task_cap_s=TASK_CAP, task_remaining_s=TASK_CAP-task_used,
            global_charge_s=charge, global_remaining_s=budget.CAP-charge,
            reporting_reserve_s=budget.RESERVE, reserve_preserved=budget.CAP-charge>=budget.RESERVE,
            numerical_workers=1, numerical_threads=1, gpu_s=0,
            ledger_sha256=hashlib.sha256((RUN/'ledger.jsonl').read_bytes()).hexdigest())
        if a.stage == 'verify':
            hashes = json.loads((HERE / 'preservation.json').read_text())
            summary['changed_preserved_files'] = [rel for rel, expected in hashes.items()
                if not (ROOT/rel).is_file() or hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()!=expected]
            summary['preserved_file_count'] = len(hashes)
            with (HERE / ('resource_audit_job_%d.json' % summary['through_job'])).open('x') as f:
                json.dump(summary, f, indent=2)
        print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
