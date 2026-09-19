"""Hash-verified stages. Promotion gates are executable dependencies."""
from common import *
import os,sys,subprocess,resource,shlex,fcntl
stages={'reuse':'verify_reuse.py','correctness':'correctness.py','pilot':None,'confirmation':None,'paper':'revise_paper.py','validate':'validate_delivery.py','reports':'reports.py','package':'package.py'}
def gate_block(stage,gates):
 required={'pilot':['G1','G2'],'confirmation':['G1','G2','G3']}.get(stage,[])
 return [g for g in required if gates[g]['status']!='PASS']
with (RUN/'runner.lock').open('a') as lock:
 fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 for stage in (sys.argv[1:] or list(stages)):
  gates=json.loads((DOC/'gate_decision.json').read_text())['gates'];blocked=gate_block(stage,gates)
  if blocked:
   write(RUN/(stage+'.stage.json'),dict(stage=stage,status='blocked',failed_dependencies=blocked,gate_sha256=sha(DOC/'gate_decision.json'),fits=0));print('BLOCKED',stage,blocked,flush=True);continue
  name=stages[stage]
  if name is None:raise RuntimeError('No reviewed positive-pilot protocol exists. Gate status alone cannot authorize unprepared training.')
  deps=[HERE/name,HERE/'common.py',RUN/'correctness_protocol.json',RUN/'correctness_protocol_freeze.json',DOC/'gate_decision.json']
  extra={'correctness':['finite_range.py'],'paper':['workflow.py'],'validate':['validate_minimal_tex.py'],'package':[p.name for p in HERE.glob('*.py')]}
  deps += [HERE/p for p in extra.get(stage,[])]
  deps += [DOC/p for p in ['novelty_matrix.md','method_and_proof.md','theorem_to_code.md','contribution_statement.md']]
  if stage=='reuse':deps += [PARENT/'protocol.json',PARENT/'fit_summary.csv',PARENT/'focused_membership.csv']+list((PARENT/'evaluation').glob('*.csv'))
  if stage=='paper':deps += [p for p in (ROOT/'papers/interaction_recoverability_iclr2027/v9/source').rglob('*') if p.is_file()]
  if stage in ['validate','reports','package']:deps+=list((PAPER/'source').rglob('*'));deps=[p for p in deps if p.is_file()]
  signature={str(p.relative_to(ROOT)):sha(p) for p in deps};done=RUN/(stage+'.stage.json')
  if done.exists():
   old=json.loads(done.read_text())
   if old.get('signature')==signature:
    for n,h in old['outputs'].items():assert sha(ROOT/n)==h,n
    if stage=='reuse':
     for record in json.loads((RUN/'reuse/evidence_reuse_manifest.json').read_text()):assert sha(ROOT/record['parent_path'])==record['sha256'],record['parent_path']
    print('VERIFIED',stage,flush=True);continue
   if stage in ['correctness','reuse']:
    # These read-only stages can be repeated after code changes; previous receipt preserved.
    (RUN/'logs'/f'{stage}-superseded-{time.time_ns()}.json').write_text(done.read_text())
   else:(RUN/'logs'/f'{stage}-superseded-{time.time_ns()}.json').write_text(done.read_text())
  dirs=[RUN,DOC,PAPER];before={str(p):(p.stat().st_size,p.stat().st_mtime_ns) for d in dirs for p in d.rglob('*') if p.is_file()}
  t=time.monotonic();cpu=resource.getrusage(resource.RUSAGE_CHILDREN);log=RUN/'logs'/f'{stage}-{time.time_ns()}.log';cmd=[sys.executable,str(HERE/name)]
  env=os.environ.copy();env.update(OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',MKL_NUM_THREADS='2')
  with log.open('w') as f:r=subprocess.run(cmd,cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT)
  c=resource.getrusage(resource.RUSAGE_CHILDREN);record=dict(stage=stage,status='executed' if r.returncode==0 else 'failed',command=shlex.join(cmd),exit_code=r.returncode,wall_s=time.monotonic()-t,child_cpu_s=c.ru_utime+c.ru_stime-cpu.ru_utime-cpu.ru_stime,max_rss_kib=c.ru_maxrss,log=str(log.relative_to(ROOT)),signature=signature)
  with (RUN/'commands.jsonl').open('a') as f:f.write(json.dumps(record)+'\n')
  print(stage,record['status'],round(record['wall_s'],3),flush=True)
  if r.returncode:print(log.read_text()[-6000:]);sys.exit(r.returncode)
  record['outputs']={str(p.relative_to(ROOT)):sha(p) for d in dirs for p in d.rglob('*') if p.is_file() and p.name not in ['commands.jsonl','runner.lock'] and before.get(str(p))!=(p.stat().st_size,p.stat().st_mtime_ns)}
  write(done,record)
