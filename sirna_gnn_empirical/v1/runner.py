"""Resumable campaign launcher. Fresh task accounting; no historical allowance use."""
import argparse,datetime,hashlib,json,os,resource,signal,subprocess,sys,time,fcntl
from pathlib import Path
CODE=Path(__file__).resolve().parent; REPO=CODE.parent.parent
os.environ['PYTHONPATH']=str(CODE/'_deps')+os.pathsep+os.environ.get('PYTHONPATH','')
STAGES=['acquire','dependencies','inspect_sources','supplements','patent_sources','primary_tables','reconcile','verify_assays','split','graphs','profile','develop','freeze','train','predict','evaluate','figures','paper']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--stage',choices=STAGES);ap.add_argument('--through',choices=STAGES);ap.add_argument('--run',default=(CODE/'active_run.txt').read_text().strip());a=ap.parse_args()
 root=(REPO/a.run).resolve();root.mkdir(parents=True,exist_ok=True)
 os.sched_setaffinity(0,set(sorted(os.sched_getaffinity(0))[:2]))
 lock=(root/'campaign.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 stages=[a.stage] if a.stage else STAGES[:STAGES.index(a.through)+1] if a.through else STAGES
 for stage in stages:
  done=root/f'{stage}.done.json'
  if done.exists():
   rec=json.loads(done.read_text())
   for item in rec.get('output_hashes',[]):
    assert sha(root/item['path'])==item['sha256'],item
   print('RESUME verified completed stage:',stage,flush=True);continue
  if (root/'resource_plan.json').exists() and (root/'commands.jsonl').exists():
   prior=[json.loads(x) for x in (root/'commands.jsonl').read_text().splitlines()];plan=json.loads((root/'resource_plan.json').read_text())
   assert sum(x['child_cpu_s'] for x in prior)<plan['maximum_campaign_child_cpu_s'],'Campaign CPU stop condition reached'
   assert sum(x['wall_s'] for x in prior)<plan['maximum_campaign_aggregate_command_wall_s'],'Campaign aggregate wall stop condition reached'
  stagefile=CODE/('training_profile.py' if stage=='profile' else stage+'.py');assert stagefile.exists(),f'Stage implementation missing: {stagefile}'
  pending=root/f'{stage}.running.json';assert not pending.exists(),f'Unclosed stage requires explicit reconciliation: {pending}'
  timestamp=datetime.datetime.now(datetime.timezone.utc).isoformat();before=resource.getrusage(resource.RUSAGE_CHILDREN);start=time.monotonic()
  command=[sys.executable,str(stagefile),'--run',str(root)]
  caps={'acquire':300,'dependencies':180,'inspect_sources':120,'supplements':120,'patent_sources':180,'primary_tables':180,'reconcile':300,'verify_assays':180,'split':120,'graphs':120,'profile':120,'develop':900,'freeze':30,'train':1800,'predict':180,'evaluate':180,'figures':120,'paper':180}
  # These are conservative task-specific stop conditions, not a grant or a reused allowance.
  record={'stage':stage,'command':command,'started_utc':timestamp,'timeout_s':caps[stage],'script_sha256':sha(stagefile),'cpu_affinity':sorted(os.sched_getaffinity(0)),'environment':{k:os.environ.get(k) for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','CUDA_VISIBLE_DEVICES']}}
  pending.write_text(json.dumps(record,indent=2)+'\n')
  log=root/'logs'/f'{stage}-{time.time_ns()}.log';log.parent.mkdir(exist_ok=True)
  with log.open('w') as out:
   proc=subprocess.Popen(command,cwd=REPO,stdout=out,stderr=subprocess.STDOUT,start_new_session=True)
   try:code=proc.wait(timeout=caps[stage]);timed_out=False
   except subprocess.TimeoutExpired:
    timed_out=True;os.killpg(proc.pid,signal.SIGTERM)
    try:code=proc.wait(timeout=10)
    except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);code=proc.wait()
  after=resource.getrusage(resource.RUSAGE_CHILDREN)
  record.update(ended_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),exit_code=code,timed_out=timed_out,wall_s=time.monotonic()-start,child_cpu_s=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,child_peak_rss_kib=after.ru_maxrss,log=str(log.relative_to(root)),measurement_scope='Waited child and descendants; peak RSS is process-family lifetime high-water mark, not additive. GPU device/synchronized active wall recorded by fitting stages; does not claim GPU kernel time.')
  outmanifest=root/f'{stage}.outputs.json'
  if outmanifest.exists():record['output_hashes']=[{'path':p,'sha256':sha(root/p)} for p in json.loads(outmanifest.read_text())]
  with (root/'commands.jsonl').open('a') as out:out.write(json.dumps(record)+'\n')
  pending.unlink()
  if code:
   (root/f'{stage}.failed-{time.time_ns()}.json').write_text(json.dumps(record,indent=2)+'\n');print(log.read_text()[-6000:]);raise SystemExit(code)
  done.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record),flush=True)
if __name__=='__main__':main()
