"""Resumable content-hashed follow-up stages; stdout and resource accounting retained."""
from common import *
import sys,subprocess,fcntl,resource,shlex
stage_map={'adjudicate':'adjudicate.py','plan':'plan.py','fit':'train_focused.py','evaluate':'evaluate.py','visibility':'visibility.py','review_checks':'review_checks.py','ranking':'ranking.py','source_checks':'source_checks.py','audit':'scientific_audit.py','paper':'paper.py','validate':'validate_delivery.py','package':'package.py'}
selected=sys.argv[1:] or list(stage_map)
with (RUN/'campaign.lock').open('a') as lock:
 fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 for stage in selected:
  script=HERE/stage_map[stage];assert script.is_file(),script
  deps=[script,HERE/'common.py']+([HERE/'engine.py',HERE/'architecture.py',RUN/'protocol.json',RUN/'protocol_freeze.json',RUN/'adjudication/primary_reanchored_targets.csv'] if stage in ['fit','evaluate','visibility'] else [])
  if stage=='paper':deps += [HERE/'appendix_author.py',HERE/'authoring/main_body.tex',HERE/'authoring/robustness_appendix.tex',RUN/'evaluation/activity_metrics.csv',RUN/'visibility/input_counts.csv']
  if stage=='validate':deps += list((ROOT/'papers/interaction_recoverability_iclr2027/v9/source').rglob('*'));deps=[p for p in deps if p.is_file()]
  signature={str(p.relative_to(ROOT)):sha(p) for p in deps};done=RUN/f'{stage}.done.json'
  if done.exists():
   rec=json.loads(done.read_text());assert rec['signature']==signature
   for n,h in rec['outputs'].items():assert sha(RUN/n)==h,n
   print('VERIFIED',stage,flush=True);continue
  before={str(p.relative_to(RUN)):(p.stat().st_size,p.stat().st_mtime_ns) for p in RUN.rglob('*') if p.is_file()};logs=RUN/'logs';logs.mkdir(exist_ok=True);log=logs/f'{stage}-{time.time_ns()}.log';t=time.monotonic();cpu=resource.getrusage(resource.RUSAGE_CHILDREN);cmd=[sys.executable,str(script)]
  with log.open('w') as out:result=subprocess.run(cmd,cwd=ROOT,stdout=out,stderr=subprocess.STDOUT)
  after=resource.getrusage(resource.RUSAGE_CHILDREN);rec=dict(stage=stage,command=shlex.join(cmd),exit_code=result.returncode,wall_s=time.monotonic()-t,child_cpu_s=after.ru_utime+after.ru_stime-cpu.ru_utime-cpu.ru_stime,max_rss_kib=after.ru_maxrss,log=str(log.relative_to(RUN)),signature=signature)
  with (RUN/'commands.jsonl').open('a') as f:f.write(json.dumps(rec)+'\n')
  print(json.dumps(rec),flush=True)
  if result.returncode:print(log.read_text()[-5000:]);sys.exit(result.returncode)
  rec['outputs']={str(p.relative_to(RUN)):sha(p) for p in RUN.rglob('*') if p.is_file() and p.name not in ['commands.jsonl','campaign.lock'] and (str(p.relative_to(RUN)) not in before or before[str(p.relative_to(RUN))]!=(p.stat().st_size,p.stat().st_mtime_ns))};write(done,rec)
