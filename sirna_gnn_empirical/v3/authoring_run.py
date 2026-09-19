"""Meter manuscript-generation attempts without advertising a frozen completed stage."""
from common import *
import sys,subprocess,resource,shlex
for name in sys.argv[1:]:
 assert name in ['figures','paper','build_audit'],name
 cmd=[sys.executable,str(HERE/(name+'.py'))];out=RUN/'authoring_attempts';out.mkdir(exist_ok=True);stamp=str(time.time_ns());t=time.perf_counter();c=resource.getrusage(resource.RUSAGE_CHILDREN);env=dict(os.environ,OMP_NUM_THREADS='2',MKL_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',PYTHONDONTWRITEBYTECODE='1')
 with (out/(name+'-'+stamp+'.log')).open('w') as f:status=subprocess.call(cmd,cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT)
 cc=resource.getrusage(resource.RUSAGE_CHILDREN);rec=dict(command=shlex.join(cmd),exit_code=status,wall_s=time.perf_counter()-t,child_cpu_s=cc.ru_utime+cc.ru_stime-c.ru_utime-c.ru_stime,code_sha256=sha(HERE/(name+'.py')),log=str((out/(name+'-'+stamp+'.log')).relative_to(RUN)),scope='Figure/table/document authoring and checking; no model fitting. Later registered stages and nested build records may repeat or overlap this administration.')
 write(out/(name+'-'+stamp+'.json'),rec);print(json.dumps(rec),flush=True)
 if status:print((out/(name+'-'+stamp+'.log')).read_text()[-8000:]);sys.exit(status)
