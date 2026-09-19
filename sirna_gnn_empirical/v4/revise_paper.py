"""Preserve each typesetting attempt and its receipt; regenerate only paper artifacts."""
from common import *
import subprocess,sys,shutil,resource,shlex
P=ROOT/'papers/interaction_recoverability_iclr2027/v9';attempt=RUN/'authoring_attempts'/str(time.time_ns());attempt.mkdir(parents=True)
for n in ['main.tex','appendix.tex']:shutil.copy2(P/'source'/n,attempt/n)
if (RUN/'paper.done.json').exists():shutil.move(RUN/'paper.done.json',attempt/'paper.done.json')
for name,cmd in [('generate',['bash',str(HERE/'run_all.sh'),'paper']),('build',['bash',str(P/'typesetting/build.sh')])]:
 t=time.monotonic();c=resource.getrusage(resource.RUSAGE_CHILDREN);result=subprocess.run(cmd,cwd=ROOT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);z=resource.getrusage(resource.RUSAGE_CHILDREN);log=attempt/(name+'.log');log.write_text(result.stdout);record=dict(stage='authoring_'+name,command=shlex.join(cmd),exit_code=result.returncode,wall_s=time.monotonic()-t,child_cpu_s=z.ru_utime+z.ru_stime-c.ru_utime-c.ru_stime,max_rss_kib=z.ru_maxrss,log=str(log.relative_to(RUN)))
 with (RUN/'authoring_commands.jsonl').open('a') as f:f.write(json.dumps(record)+'\n')
 print(result.stdout[-4000:],flush=True)
 if result.returncode:raise SystemExit(result.returncode)
