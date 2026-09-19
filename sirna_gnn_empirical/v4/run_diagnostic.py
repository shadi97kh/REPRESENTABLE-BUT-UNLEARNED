"""Record a disjoint read-only-input diagnostic while the sole fitter runs."""
from common import *
import subprocess,sys,resource,shlex
stage=sys.argv[1];assert stage in ['review_checks','ranking','source_checks'];script=HERE/f'{stage}.py';sig={str(p.relative_to(ROOT)):sha(p) for p in [script,HERE/'common.py']};done=RUN/f'{stage}.done.json'
if done.exists():
 rec=json.loads(done.read_text());assert rec['signature']==sig
 for n,h in rec['outputs'].items():assert sha(RUN/n)==h
 print('VERIFIED',stage);raise SystemExit
outdir=RUN/stage;outdir.mkdir(exist_ok=True);log=RUN/'logs'/f'{stage}-{time.time_ns()}.log';t=time.monotonic();c=resource.getrusage(resource.RUSAGE_CHILDREN);cmd=[sys.executable,str(script)]
with log.open('w') as f:result=subprocess.run(cmd,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
z=resource.getrusage(resource.RUSAGE_CHILDREN);rec=dict(stage=stage,command=shlex.join(cmd),exit_code=result.returncode,wall_s=time.monotonic()-t,child_cpu_s=z.ru_utime+z.ru_stime-c.ru_utime-c.ru_stime,max_rss_kib=z.ru_maxrss,log=str(log.relative_to(RUN)),signature=sig)
with (RUN/'diagnostic_commands.jsonl').open('a') as f:f.write(json.dumps(rec)+'\n')
if result.returncode:print(log.read_text()[-3000:]);raise SystemExit(result.returncode)
rec['outputs']={str(p.relative_to(RUN)):sha(p) for p in outdir.rglob('*') if p.is_file()};write(done,rec);print(log.read_text()[-3000:])
