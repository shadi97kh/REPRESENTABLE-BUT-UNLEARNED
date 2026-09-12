"""Serial one-core process-group watchdog and append-only preparation ledger."""
import os
import time
START = time.monotonic()
CPU = min(os.sched_getaffinity(0))
os.sched_setaffinity(0, {CPU})
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
os.environ['CUDA_VISIBLE_DEVICES']=''
import argparse
import fcntl
import json
from pathlib import Path
import resource
import signal
import subprocess
import sys


def append(path, event):
    with path.open('a') as f:
        f.write(json.dumps(event,allow_nan=False)+'\n'); f.flush(); os.fsync(f.fileno())


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-root',required=True)
    p.add_argument('--init',action='store_true')
    p.add_argument('command',nargs=argparse.REMAINDER)
    a=p.parse_args(); root=Path(a.run_root)
    if a.init:
        root.mkdir(parents=True,exist_ok=False)
        append(root/'ledger.jsonl',dict(event='allowance',cpu_s_cap=600,aggregate_job_wall_s_cap=600,
            gpu_s_cap=0,workers=1,numerical_threads=1,affinity_cpu=CPU,
            reserve_s=30,conservative_charge_s=30,measured_cpu_s=None,
            note='Initial text-only repository/prompt discovery preceded this meter. Thirty seconds reserved/charged for that unmeasured setup and launcher overhead; not an instrumented CPU measurement. All subsequent local jobs use inherited one-core affinity. Historical balance unknown.',
            started_utc='2026-09-12T08:21:54Z'))
        print(root);return
    lock=(root/'lock').open('a')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    events=[json.loads(x) for x in (root/'ledger.jsonl').read_text().splitlines()]
    if any(e['event']=='start' and not any(t.get('job_id')==e['job_id'] and t['event']=='end' for t in events) for e in events):
        raise RuntimeError('Unclosed job: reconcile before continuing')
    used=sum(e.get('conservative_charge_s',0) for e in events)
    remaining=600-used-30  # retain a clean-stop reserve
    if remaining<=0:raise RuntimeError('Preparation allowance exhausted')
    command=a.command[1:] if a.command[:1]==['--'] else a.command
    if not command:raise ValueError('Explicit command required')
    job_id=sum(e['event']=='start' for e in events)+1
    append(root/'ledger.jsonl',dict(event='start',job_id=job_id,command=command,timeout_s=remaining,
           affinity=[CPU],environment={k:os.environ[k] for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','CUDA_VISIBLE_DEVICES')}))
    before=resource.getrusage(resource.RUSAGE_CHILDREN)
    proc=subprocess.Popen(command,start_new_session=True)
    timed_out=False
    try:code=proc.wait(timeout=remaining)
    except subprocess.TimeoutExpired:
        timed_out=True;os.killpg(proc.pid,signal.SIGKILL);code=proc.wait()
    finally:
        # No descendant is allowed to persist beyond the job, even on success.
        try:os.killpg(proc.pid,signal.SIGKILL)
        except ProcessLookupError:pass
    child=resource.getrusage(resource.RUSAGE_CHILDREN)
    parent=resource.getrusage(resource.RUSAGE_SELF)
    wall=time.monotonic()-START
    measured=parent.ru_utime+parent.ru_stime+child.ru_utime+child.ru_stime-before.ru_utime-before.ru_stime
    charge=max(wall,measured)+.05
    append(root/'ledger.jsonl',dict(event='end',job_id=job_id,returncode=code,timed_out=timed_out,
        wall_s=wall,measured_cpu_s=measured,conservative_charge_s=charge,
        bound='Whole job restricted to one CPU shared with descendants; elapsed wall plus launcher margin bounds concurrent tree CPU. Escaping affinity/process group is prohibited.',
        measurement_limit='Waited-child CPU plus parent CPU; un-reaped grandchildren may be absent from measured counters but covered by one-core wall bound.',
        parent_peak_rss_kib=parent.ru_maxrss,waited_child_peak_rss_kib=child.ru_maxrss))
    final_events=[json.loads(line) for line in (root/'ledger.jsonl').read_text().splitlines()]
    ends=[e for e in final_events if e['event']=='end']
    summary=dict(through_job=job_id,measured_cpu_s=sum(e['measured_cpu_s'] for e in ends),
        measured_cpu_core_minutes=sum(e['measured_cpu_s'] for e in ends)/60,
        aggregate_instrumented_job_wall_s=sum(e['wall_s'] for e in ends),
        conservative_charge_s=sum(e.get('conservative_charge_s',0) for e in final_events),
        initial_unmeasured_setup_charge_s=30,per_job_launcher_margin_s=.05,
        cpu_s_cap=600,aggregate_job_wall_s_cap=600,gpu_s=0,workers=1,numerical_threads=1,
        scope='Fresh preparation allowance only; initial text inspection charged but unmetered; one-core wall bound covers all subsequent local jobs and descendants',
        task_started_utc='2026-09-12T08:21:54Z',
        budget_remaining_s=600-used-charge,
        peak_waited_child_rss_kib=max(e['waited_child_peak_rss_kib'] for e in ends),
        peak_parent_rss_kib=max(e['parent_peak_rss_kib'] for e in ends))
    from datetime import datetime,timezone
    summary['recorded_utc']=datetime.now(timezone.utc).isoformat()
    summary['elapsed_task_wall_s']=(datetime.now(timezone.utc)-datetime.fromisoformat('2026-09-12T08:21:54+00:00')).total_seconds()
    with (root/f'resource_snapshot-{job_id}.json').open('x') as out:json.dump(summary,out,indent=2)
    print(json.dumps(summary))
    raise SystemExit(code)


if __name__=='__main__':main()
