"""Wait for the active stage, then execute the remaining frozen scientific stages."""
from common import *
import fcntl,subprocess,sys
with (RUN/'campaign.lock').open('a') as f:
 fcntl.flock(f,fcntl.LOCK_EX)
 assert (RUN/'activity.done.json').exists(),'Activity did not finish successfully; no downstream stage is launched.'
 fcntl.flock(f,fcntl.LOCK_UN)
for stage in ['pair','b3_evaluate','evaluate','analysis_supplement']:
 print('EXECUTING',stage,flush=True)
 subprocess.run(['bash',str(HERE/'run_all.sh'),stage],cwd=ROOT,check=True)
