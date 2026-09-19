"""Local attainment task cap; reuses unchanged current-ledger single-core watchdog."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).resolve().parent
RUN=ROOT/'runs/interaction_recoverability/learned-target-v1-20260912-164000'
TASK='local-attainment-v1'
EXPECTED='7f30adc8bec106fda1801663a57a31c3578578a2d0259a073a27639f65925225'
sys.path.insert(0,str(ROOT));sys.dont_write_bytecode=True
os.environ['PYTHONDONTWRITEBYTECODE']='1'
from interaction_recoverability.learned_target import budget

def main():
    p=argparse.ArgumentParser();p.add_argument('--seconds',type=float,default=45.)
    p.add_argument('stage',choices=['register','run','verify']);p.add_argument('command',nargs=argparse.REMAINDER)
    a=p.parse_args();lock=(HERE/'task.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    rows=budget.events(RUN)
    starts={r['job_id'] for r in rows if r['event']=='start'}
    ends={r['job_id'] for r in rows if r['event']=='end'}
    assert starts==ends,'Unclosed job'
    assert rows[0]['cpu_s_cap']==rows[0]['aggregate_job_wall_s_cap']==7200
    assert rows[0]['reserve_s']==600
    config=json.loads((RUN/'frozen_config.json').read_text())
    for rel,sha in config['learner_code_sha256'].items():
        assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==sha,rel
    used=sum(r.get('conservative_charge_s',0.) for r in rows)
    markers=[r for r in rows if r.get('event')=='continuation' and r.get('task')==TASK]
    if a.stage=='register':
        assert not markers,'No reset permitted'
        assert hashlib.sha256((RUN/'ledger.jsonl').read_bytes()).hexdigest()==EXPECTED
        assert 7200-used-600>=600
        paths=[q for q in ROOT.iterdir() if q.is_file()]
        for base in ['docs/interaction_recoverability','interaction_recoverability','configs','tests',str(RUN.relative_to(ROOT))]:
            paths.extend(q for q in (ROOT/base).rglob('*') if q.is_file() and '__pycache__' not in q.parts)
        hashes={str(q.relative_to(ROOT)):hashlib.sha256(q.read_bytes()).hexdigest() for q in paths
                if HERE not in q.parents and q.name not in ['ledger.jsonl','lock','task.lock']}
        old=ROOT/'runs/interaction_recoverability/prep-20260912-082154/ledger.jsonl'
        hashes[str(old.relative_to(ROOT))]=hashlib.sha256(old.read_bytes()).hexdigest()
        with (HERE/'preservation.json').open('x') as f:json.dump(hashes,f,indent=2)
        budget.append(RUN/'ledger.jsonl',dict(event='continuation',task=TASK,
            additional_cpu_s_cap=600.,additional_job_wall_s_cap=600.,prior_charge_s=used,
            prior_through_job=max(ends),prior_ledger_sha256=EXPECTED,conservative_charge_s=150.,
            measured_cpu_s=None,gpu_s_cap=0,workers=1,numerical_threads=1,
            note='Bounded local attainment continuation within the existing allowance; no new allowance. 150s conservative debit for bootstrap, static research/agent reads, file authoring and finalization outside jobs. No allowance reset or historical-ledger write.'))
        rows=budget.events(RUN);markers=[r for r in rows if r.get('event')=='continuation' and r.get('task')==TASK]
    assert len(markers)==1,'Registration required'
    marker=markers[0];used=sum(r.get('conservative_charge_s',0.) for r in rows)
    task_used=used-marker['prior_charge_s'];remaining=min(600-task_used,7200-used-600)
    if a.stage=='run':
        assert 0<a.seconds<=120
        timeout=min(a.seconds,remaining-5.)
        assert timeout>0,'No admitted task allowance'
        cmd=a.command[1:] if a.command[:1]==['--'] else a.command
        assert cmd
        sys.argv=['unchanged-budget','--run-root',str(RUN),'--phase','compute','--job-seconds',str(timeout),'--',*cmd]
        budget.main()
    else:
        completed=[r for r in rows if r['event']=='end' and r['job_id']>marker['prior_through_job']]
        summary=dict(task=TASK,through_job=max(ends),task_cpu_s=sum(r['measured_cpu_s'] for r in completed),
            task_wall_s=sum(r['wall_s'] for r in completed),task_charge_s=task_used,task_cap_s=600.,
            task_remaining_s=600-task_used,global_charge_s=used,global_remaining_s=7200-used,
            reporting_reserve_s=600.,reserve_preserved=7200-used>=600.,gpu_s=0,workers=1,threads=1,
            ledger_sha256=hashlib.sha256((RUN/'ledger.jsonl').read_bytes()).hexdigest())
        if a.stage=='verify':
            hashes=json.loads((HERE/'preservation.json').read_text())
            summary['preserved_file_count']=len(hashes)
            summary['changed_preserved_files']=[rel for rel,sha in hashes.items() if not (ROOT/rel).is_file() or hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()!=sha]
            summary['job_statuses']=[dict(job_id=r['job_id'],returncode=r['returncode'],timed_out=r['timed_out'],affinity_violations=r['affinity_violations']) for r in completed]
            with (HERE/('resource_audit_job_%d.json'%max(ends))).open('x') as f:json.dump(summary,f,indent=2)
        print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
