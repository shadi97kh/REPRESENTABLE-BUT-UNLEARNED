"""Nested 120s efficiency cap within existing practical 300s and global reserve."""
import argparse, fcntl, hashlib, json, os, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).resolve().parent
RUN=ROOT/'runs/interaction_recoverability/learned-target-v1-20260912-164000'
TASK='efficiency-audit-v1'
PARENT='practical-feasibility-v1'
EXPECTED='f446782bef0bf20ac16c99bf801c2f1a5e0d569959cb01c54e63982a1ceb6a96'
sys.dont_write_bytecode=True
os.environ['PYTHONDONTWRITEBYTECODE']='1'
sys.path.insert(0,str(ROOT))
from interaction_recoverability.learned_target import budget

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    p=argparse.ArgumentParser();p.add_argument('--seconds',type=float,default=30)
    p.add_argument('stage',choices=['register','run','verify'])
    p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args()
    lock=(HERE/'task.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    rows=budget.events(RUN)
    starts={r['job_id'] for r in rows if r['event']=='start'}
    ends={r['job_id'] for r in rows if r['event']=='end'}
    assert starts==ends,'Unclosed job'
    assert rows[0]['cpu_s_cap']==rows[0]['aggregate_job_wall_s_cap']==7200
    assert rows[0]['reserve_s']==600
    frozen=json.loads((RUN/'frozen_config.json').read_text())
    for rel,h in frozen['learner_code_sha256'].items():assert sha(ROOT/rel)==h,rel
    parent=[r for r in rows if r.get('task')==PARENT and r['event']=='continuation']
    assert len(parent)==1 and parent[0]['additional_cpu_s_cap']==parent[0]['additional_job_wall_s_cap']==300
    parent=parent[0]
    used=sum(r.get('conservative_charge_s',0) for r in rows)
    markers=[r for r in rows if r.get('task')==TASK and r['event']=='continuation']
    if a.stage=='register':
        assert not markers,'No reset allowed'
        assert sha(RUN/'ledger.jsonl')==EXPECTED
        assert min(300-(used-parent['prior_charge_s']),7200-used-600)>=120
        paths=[q for q in ROOT.iterdir() if q.is_file()]
        for base in ['docs','paper','interaction_recoverability','intervention','scmp','certmp','experiments','configs','tests','runs']:
            paths.extend(q for q in (ROOT/base).rglob('*') if q.is_file() and '__pycache__' not in q.parts)
        hashes={str(q.relative_to(ROOT)):sha(q) for q in paths
                if HERE not in q.parents and q.name not in ['ledger.jsonl','lock','task.lock']}
        old=ROOT/'runs/interaction_recoverability/prep-20260912-082154/ledger.jsonl'
        hashes[str(old.relative_to(ROOT))]=sha(old)
        with (HERE/'preservation.json').open('x') as f:json.dump(hashes,f,indent=2)
        budget.append(RUN/'ledger.jsonl',dict(event='continuation',task=TASK,parent_task=PARENT,
            additional_cpu_s_cap=120.,additional_job_wall_s_cap=120.,prior_charge_s=used,
            prior_through_job=max(ends),prior_ledger_sha256=EXPECTED,conservative_charge_s=60.,
            measured_cpu_s=None,gpu_s_cap=0,workers=1,numerical_threads=1,
            frozen_calculation_sha256=sha(HERE/'calculate.py'),
            note='Nested audit; also charged to practical-feasibility 300s. 60s conservative debit for bootstrap, static reads/internal reviews, authoring and finalization outside jobs. No allowance reset. Twelve directions and two quadrature configurations frozen before outputs.'))
        rows=budget.events(RUN);markers=[r for r in rows if r.get('task')==TASK and r['event']=='continuation']
    assert len(markers)==1
    marker=markers[0]
    assert sha(HERE/'calculate.py')==marker['frozen_calculation_sha256'],'Frozen calculation modified'
    used=sum(r.get('conservative_charge_s',0) for r in rows)
    task_used=used-marker['prior_charge_s'];parent_used=used-parent['prior_charge_s']
    remaining=min(120-task_used,300-parent_used,7200-used-600)
    if a.stage=='run':
        assert 0<a.seconds<=30
        timeout=min(a.seconds,remaining-5.)
        assert timeout>0
        cmd=a.command[1:] if a.command[:1]==['--'] else a.command
        assert cmd
        sys.argv=['unchanged-budget','--run-root',str(RUN),'--phase','compute','--job-seconds',str(timeout),'--',*cmd]
        budget.main()
    else:
        completed=[r for r in rows if r['event']=='end' and r['job_id']>marker['prior_through_job']]
        summary=dict(task=TASK,through_job=max(ends),measured_cpu_s=sum(r['measured_cpu_s'] for r in completed),
            aggregate_job_wall_s=sum(r['wall_s'] for r in completed),overhead_debit_s=60.,
            conservative_charge_s=task_used,task_cap_s=120.,task_remaining_s=120-task_used,
            parent_charge_s=parent_used,parent_cap_s=300.,parent_remaining_s=300-parent_used,
            global_charge_s=used,global_remaining_s=7200-used,reporting_reserve_s=600.,
            reserve_preserved=7200-used>=600,gpu_s=0,workers=1,threads=1,
            ledger_sha256=sha(RUN/'ledger.jsonl'))
        if a.stage=='verify':
            hashes=json.loads((HERE/'preservation.json').read_text())
            summary['preserved_file_count']=len(hashes)
            summary['changed_preserved_files']=[rel for rel,h in hashes.items() if not (ROOT/rel).is_file() or sha(ROOT/rel)!=h]
            summary['job_records']=[r for r in rows if r.get('job_id',0)>marker['prior_through_job']]
            summary['artifact_hashes']={q.name:sha(q) for q in HERE.iterdir() if q.is_file() and q.name not in ['resources.json','task.lock','commands_and_resources.md']}
            assert not summary['changed_preserved_files']
            assert task_used<=120 and parent_used<=300 and summary['reserve_preserved']
            with (HERE/'resources.json').open('x') as f:json.dump(summary,f,indent=2)
        print(json.dumps(summary,indent=2))
if __name__=='__main__':main()

