"""New 7200-second allowance, one-core tree watchdog, append-only accounting."""
import os
import time
START=time.monotonic()
CPU=min(os.sched_getaffinity(0))
os.sched_setaffinity(0,{CPU})
for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
    os.environ[name]='1'
os.environ['CUDA_VISIBLE_DEVICES']=''
import argparse
import fcntl
import hashlib
import json
from pathlib import Path
import resource
import signal
import subprocess
import sys
from datetime import datetime,timezone

CAP=7200.
RESERVE=600.

def append(path,event):
    event={'utc':datetime.now(timezone.utc).isoformat(),**event}
    with path.open('a') as f:
        f.write(json.dumps(event,allow_nan=False)+'\n');f.flush();os.fsync(f.fileno())

def events(root):
    return [json.loads(x) for x in (root/'ledger.jsonl').read_text().splitlines()]

def tree_pids(parent):
    """Observe the process tree, including children which changed process groups."""
    parents={}
    for p in Path('/proc').iterdir():
        if p.name.isdigit():
            try:
                raw=(p/'stat').read_text();fields=raw[raw.rfind(')')+2:].split()
                parents[int(p.name)]=int(fields[1])
            except (OSError,ValueError,IndexError):pass
    found={parent}
    while True:
        new={pid for pid,ppid in parents.items() if ppid in found}
        if new<=found:return found
        found|=new

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-root',required=True)
    p.add_argument('--init',action='store_true')
    p.add_argument('--phase',choices=['compute','report'],default='compute')
    p.add_argument('--job-seconds',type=float,default=600.)
    p.add_argument('command',nargs=argparse.REMAINDER)
    a=p.parse_args();root=Path(a.run_root)
    if a.init:
        root.mkdir(parents=True,exist_ok=False)
        old=Path('runs/interaction_recoverability/prep-20260912-082154/ledger.jsonl')
        append(root/'ledger.jsonl',dict(event='allowance',cpu_s_cap=CAP,aggregate_job_wall_s_cap=CAP,
            reserve_s=RESERVE,gpu_s_cap=0,workers=1,numerical_threads=1,conservative_charge_s=30.,
            measured_cpu_s=None,historical_ledger=str(old),historical_sha256=hashlib.sha256(old.read_bytes()).hexdigest(),
            note='Fresh explicit authorization in Codex_Learned_Interaction_Contribution_Execution.md. Thirty seconds charged for text discovery/bootstrap, unmeasured. First discovery was mistakenly logged as historical job 27 (0.101903 s charge); retained transparently and also covered here. No further historical debit.'))
        print(str(root));return
    lock=(root/'lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    prior=events(root)
    ends={e['job_id'] for e in prior if e['event']=='end'}
    if any(e['event']=='start' and e['job_id'] not in ends for e in prior):
        raise RuntimeError('Unclosed job; reconcile before admission')
    used=sum(e.get('conservative_charge_s',0) for e in prior)
    allowance=CAP-used-(RESERVE if a.phase=='compute' else 2.)
    timeout=min(a.job_seconds,allowance-2.)
    if not 0<timeout or not 0<a.job_seconds<=3600:raise RuntimeError('No admissible job allowance')
    command=a.command[1:] if a.command[:1]==['--'] else a.command
    if not command:raise ValueError('Explicit command required')
    if a.phase=='report' and not any(x in command for x in ['evaluate','report','verify']):
        raise ValueError('Reporting reserve only admits named evaluate/report/verify stages')
    job=sum(e['event']=='start' for e in prior)+1
    append(root/'ledger.jsonl',dict(event='start',job_id=job,phase=a.phase,command=command,timeout_s=timeout,affinity=[CPU]))
    before=resource.getrusage(resource.RUSAGE_CHILDREN)
    proc=subprocess.Popen(command,start_new_session=True)
    observed={proc.pid};timed_out=False;affinity_violations=[]
    try:
        while proc.poll() is None:
            observed|=tree_pids(proc.pid)
            for pid in observed:
                try:
                    if os.sched_getaffinity(pid)!={CPU}:
                        affinity_violations.append(pid);os.sched_setaffinity(pid,{CPU})
                except ProcessLookupError:pass
            if time.monotonic()-START>=timeout:
                timed_out=True;break
            time.sleep(.1)
    finally:
        # Gentle stop permits checkpoint flush; hard stop prevents detached descendants.
        if proc.poll() is None:
            try:os.killpg(proc.pid,signal.SIGTERM)
            except ProcessLookupError:pass
            try:proc.wait(timeout=.5)
            except subprocess.TimeoutExpired:pass
        observed|=tree_pids(proc.pid)
        for pid in observed:
            try:os.kill(pid,signal.SIGKILL)
            except ProcessLookupError:pass
        try:os.killpg(proc.pid,signal.SIGKILL)
        except ProcessLookupError:pass
        code=proc.wait()
    child=resource.getrusage(resource.RUSAGE_CHILDREN);parent=resource.getrusage(resource.RUSAGE_SELF)
    cpu=parent.ru_utime+parent.ru_stime+child.ru_utime+child.ru_stime-before.ru_utime-before.ru_stime
    wall=time.monotonic()-START;charge=max(cpu,wall)+.1
    append(root/'ledger.jsonl',dict(event='end',job_id=job,returncode=code,timed_out=timed_out,
        measured_cpu_s=cpu,wall_s=wall,conservative_charge_s=charge,observed_process_count=len(observed),
        affinity_violations=sorted(set(affinity_violations)),parent_peak_rss_kib=parent.ru_maxrss,
        waited_child_peak_rss_kib=child.ru_maxrss,
        tree_bound='Inherited single-core affinity makes aggregate wall a conservative entire-tree CPU bound. Waited CPU can omit unreaped grandchildren. Affinity escape is prohibited and checked; all observed descendants killed on exit.'))
    current=events(root);completed=[e for e in current if e['event']=='end']
    summary=dict(through_job=job,measured_cpu_s=sum(e['measured_cpu_s'] for e in completed),
        aggregate_job_wall_s=sum(e['wall_s'] for e in completed),
        conservative_charge_s=sum(e.get('conservative_charge_s',0) for e in current),
        cpu_s_cap=CAP,aggregate_job_wall_s_cap=CAP,reserve_s=RESERVE,gpu_s=0,workers=1,numerical_threads=1,
        recorded_utc=datetime.now(timezone.utc).isoformat())
    summary['remaining_s']=CAP-summary['conservative_charge_s']
    summary['compute_remaining_before_reserve_s']=max(0,summary['remaining_s']-RESERVE)
    with (root/f'resource_snapshot-{job}.json').open('x') as f:json.dump(summary,f,indent=2)
    print(json.dumps(summary));raise SystemExit(code if not timed_out else 124)

if __name__=='__main__':main()
