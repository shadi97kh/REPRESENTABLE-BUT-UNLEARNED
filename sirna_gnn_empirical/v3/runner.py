"""Resumable stage execution; validate code, dependency and output hashes."""
from common import *
import sys,subprocess,resource,fcntl,shlex
STAGES=['legacy_analysis','source_checks','protocol','checks','activity','pair','b3_evaluate','evaluate','analysis_supplement','parameter_audit','figures','paper','deliver']
DEPENDENCIES={'legacy_analysis':[],'source_checks':[],'protocol':['legacy_analysis'],'checks':['protocol'],'activity':['checks'],'pair':['checks'],'b3_evaluate':['activity'],'evaluate':['activity','pair','legacy_analysis','b3_evaluate'],'analysis_supplement':['evaluate'],'parameter_audit':['analysis_supplement'],'figures':['parameter_audit'],'paper':['figures','source_checks'],'deliver':['paper']}
selected=sys.argv[1:] or STAGES
lock=(RUN/'campaign.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
for stage in selected:
 if stage not in STAGES:raise ValueError(stage)
 script=HERE/f'{stage}.py';assert script.exists(),script
 deps={d:sha(RUN/f'{d}.done.json') for d in DEPENDENCIES[stage]}
 # Stage-specific source dependency list avoids invalidation from unrelated paper authoring.
 names=[f'{stage}.py','common.py']+(['metrics.py'] if stage in ['legacy_analysis','evaluate','b3_evaluate','analysis_supplement','parameter_audit'] else [])+(['engine.py','architecture.py'] if stage in ['activity','pair','checks','b3_evaluate','analysis_supplement','parameter_audit'] else [])
 codes={n:sha(HERE/n) for n in names};signature={'code':codes,'dependencies':deps};done=RUN/f'{stage}.done.json'
 if done.exists():
  old=json.loads(done.read_text());assert old['signature']==signature,(stage,'input/code change requires a new explicit stage version')
  for p,h in old['outputs'].items():assert sha(RUN/p)==h,p
  print('VERIFIED',stage,flush=True);continue
 before={str(p.relative_to(RUN)):(p.stat().st_mtime_ns,p.stat().st_size) for p in RUN.rglob('*') if p.is_file()};t=time.perf_counter();c=resource.getrusage(resource.RUSAGE_CHILDREN)
 logs=RUN/'logs';logs.mkdir(exist_ok=True);log=logs/f'{stage}-{time.time_ns()}.log';cmd=[sys.executable,str(script)];env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OMP_NUM_THREADS='2',MKL_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',CUBLAS_WORKSPACE_CONFIG=':4096:8')
 with log.open('w') as out:code=subprocess.call(cmd,cwd=ROOT,env=env,stdout=out,stderr=subprocess.STDOUT)
 cc=resource.getrusage(resource.RUSAGE_CHILDREN);rec=dict(stage=stage,command=shlex.join(cmd),wall_s=time.perf_counter()-t,child_cpu_s=cc.ru_utime+cc.ru_stime-c.ru_utime-c.ru_stime,max_child_rss_kib=cc.ru_maxrss,exit_code=code,log=str(log.relative_to(RUN)),signature=signature)
 with (RUN/'commands.jsonl').open('a') as f:f.write(json.dumps(rec)+'\n')
 print(json.dumps(rec),flush=True)
 if code:print(log.read_text()[-9000:],flush=True);sys.exit(code)
 outputs={str(p.relative_to(RUN)):sha(p) for p in RUN.rglob('*') if p.is_file() and p.name not in ['commands.jsonl','campaign.lock'] and (str(p.relative_to(RUN)) not in before or before[str(p.relative_to(RUN))]!=(p.stat().st_mtime_ns,p.stat().st_size))}
 write(done,dict(signature=signature,outputs=outputs))
