import os,sys,json,time,subprocess,resource,fcntl,hashlib
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];RUN=ROOT/(HERE/'active_run.txt').read_text().strip()
stages=['audit','source_checks','plan','factorial','grouped','b2','external','external_uncertainty','evaluate','figures','paper','integrity']
requested=sys.argv[1:] or stages
lock=(RUN/'campaign.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
for stage in requested:
 if stage not in stages:raise ValueError(stage)
 done=RUN/f'{stage}.done.json'
 if done.exists():
  d=json.loads(done.read_text())
  for p,h in d['outputs'].items():assert hashlib.sha256((RUN/p).read_bytes()).hexdigest()==h,(stage,p)
  print('VERIFIED',stage,flush=True);continue
 prior=[json.loads(s) for s in (RUN/'commands.jsonl').read_text().splitlines()] if (RUN/'commands.jsonl').exists() else []
 if stage not in ['paper','integrity']:
  assert sum(x['child_cpu_s'] for x in prior)<14400 and sum(x['wall_s'] for x in prior)<14400, 'campaign resource bound reached'
 script=HERE/(stage+'.py');cmd=[sys.executable,str(script)];env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OMP_NUM_THREADS='2',MKL_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2')
 before={str(p.relative_to(RUN)):(p.stat().st_mtime_ns,p.stat().st_size) for p in RUN.rglob('*') if p.is_file()};start=time.perf_counter();cpu=resource.getrusage(resource.RUSAGE_CHILDREN);log=RUN/'logs'/f'{stage}-{time.time_ns()}.log';log.parent.mkdir(exist_ok=True)
 with log.open('w') as f:
  proc=subprocess.Popen(cmd,cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
  try:code=proc.wait(timeout=3600)
  except subprocess.TimeoutExpired:
   import signal
   os.killpg(proc.pid,signal.SIGTERM);code=proc.wait();code=124
 aftercpu=resource.getrusage(resource.RUSAGE_CHILDREN);entry=dict(stage=stage,command=cmd,wall_s=time.perf_counter()-start,child_cpu_s=aftercpu.ru_utime+aftercpu.ru_stime-cpu.ru_utime-cpu.ru_stime,max_child_rss_kib=aftercpu.ru_maxrss,exit_code=code,log=str(log.relative_to(RUN)))
 with (RUN/'commands.jsonl').open('a') as f:f.write(json.dumps(entry)+'\n')
 print(entry,flush=True)
 if code:print(log.read_text()[-10000:]);sys.exit(code)
 outputs={str(p.relative_to(RUN)):hashlib.sha256(p.read_bytes()).hexdigest() for p in RUN.rglob('*') if p.is_file() and str(p.relative_to(RUN)) not in ['commands.jsonl','campaign.lock'] and (str(p.relative_to(RUN)) not in before or before[str(p.relative_to(RUN))]!=(p.stat().st_mtime_ns,p.stat().st_size))}
 done.write_text(json.dumps(dict(stage=stage,outputs=outputs),indent=2)+'\n')
