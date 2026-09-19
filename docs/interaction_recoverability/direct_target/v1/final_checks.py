"""Static final artifact/preservation check; no numerical-module imports."""
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
RUN=ROOT/'runs/interaction_recoverability/learned-target-v1-20260912-164000'

required=['exact_problem.md','target_identity.md','modulus.md','estimator.md','novelty.md','decision.md','commands_and_resources.md']
assert all((HERE/x).is_file() and (HERE/x).stat().st_size>100 for x in required)
preserved=json.loads((HERE/'preservation.json').read_text())
changed=[rel for rel,sha in preserved.items() if not (ROOT/rel).is_file() or hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()!=sha]
assert not changed,changed
diagnostics=json.loads((HERE/'checks.json').read_text())
assert diagnostics['status']=='all deterministic checks passed'
rows=[json.loads(x) for x in (RUN/'ledger.jsonl').read_text().splitlines()]
marker=next(r for r in rows if r.get('event')=='continuation' and r.get('task')=='direct-target-v1')
jobs=[r for r in rows if r['event']=='start' and r['job_id']>marker['prior_through_job']]
ends=[r for r in rows if r['event']=='end' and r['job_id']>marker['prior_through_job']]
assert all(r['returncode']==0 and not r['timed_out'] and not r['affinity_violations'] for r in ends)
assert sum(r.get('conservative_charge_s',0) for r in rows)-marker['prior_charge_s']<600
out={'preserved_files_checked':len(preserved),'changed_preserved_files':changed,
     'required_documents':required,'completed_task_job_ids':[r['job_id'] for r in ends],
     'active_final_check_job_ids':[r['job_id'] for r in jobs if r['job_id'] not in {e['job_id'] for e in ends}],
     'note':'Final check reads existing artifacts only. Its own completed cost is in the subsequent ledger end and snapshot.'}
with (HERE/'final_checks.json').open('x') as f:json.dump(out,f,indent=2)
with (HERE/'executed_jobs.json').open('x') as f:json.dump({'starts':jobs,'completed_ends_before_this_check':ends},f,indent=2)
print(json.dumps(out,indent=2))
