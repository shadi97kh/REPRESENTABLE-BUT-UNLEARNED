"""Bounded primary-source and official-baseline acquisition, no wet-lab work."""
from common import *
import subprocess,sys,resource,requests,shutil,importlib.util
out=RUN/'sources';out.mkdir(exist_ok=True);repo=ROOT/'handoffs/codex_revision_v1/ENsiRNA_official';deps=HERE/'_deps'
def command(label,cmd,cwd=ROOT,timeout=180):
 target=out/f'{label}.record.json'
 if target.exists():return json.loads(target.read_text())
 env=dict(os.environ,PYTHONPATH=str(deps),PYTHONDONTWRITEBYTECODE='1',OMP_NUM_THREADS='2',MKL_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2');t=time.perf_counter();c=resource.getrusage(resource.RUSAGE_CHILDREN)
 try:
  p=subprocess.run(cmd,cwd=cwd,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=timeout);code=p.returncode;log=p.stdout
 except subprocess.TimeoutExpired as e:code=124;log=(e.stdout or b'').decode() if isinstance(e.stdout,bytes) else (e.stdout or '')
 cc=resource.getrusage(resource.RUSAGE_CHILDREN);(out/f'{label}.log').write_text(log);rec=dict(command=cmd,cwd=str(cwd),exit_code=code,wall_s=time.perf_counter()-t,child_cpu_s=cc.ru_utime+cc.ru_stime-c.ru_utime-c.ru_stime,max_child_rss_kib=cc.ru_maxrss,log=f'{label}.log');write(target,rec);print(label,code,log[-1200:],flush=True);return rec
commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip();write(out/'official_pin.json',dict(repository='https://github.com/tanwenchong/ENsiRNA',commit=commit,license_statement='README Apache v2.0',code_files={str(p.relative_to(repo)):sha(p) for p in (repo/'ENsiRNA-mod').rglob('*.py')}))
command('install_dependencies',[sys.executable,'-m','pip','install','--target',str(deps),'gdown','tensorboard'])
pdb=out/'ensirna_mod_pdb.zip';command('acquire_official_pdb',[sys.executable,'-m','gdown','https://drive.google.com/uc?id=1F7cNJXMNPSjFb0UvDkDRHTkt9Tt4EGWe','-O',str(pdb),'--no-cookies'],timeout=180)
# A help invocation actually imports the official data/model dependencies; it never fits an efficacy checkpoint.
command('official_preprocess_help',[sys.executable,'-m','data.get_pdb','--help'],cwd=repo/'ENsiRNA-mod')
# Probe the exact required geometry executable without replacing coordinates.
exe='/app/rosetta/rosetta.binary.linux.release-371/main/source/bin/rna_denovo.static.linuxgccrelease'
try:command('rosetta_probe',[exe,'-help'],timeout=20)
except FileNotFoundError as e:write(out/'rosetta_probe.record.json',dict(command=[exe,'-help'],error=str(e),exit_code=None,reason='Required geometry generator is not installed at official configured path'))
acq=json.loads((out/'primary_acquisition.json').read_text()) if (out/'primary_acquisition.json').exists() else []
for name,url in [('davis_current.html','https://academic.oup.com/nar/article/53/12/gkaf479/8171869'),('davis_crossref.json','https://api.crossref.org/works/10.1093/nar/gkaf479'),('ensirna_pdb_share.html','https://drive.google.com/file/d/1F7cNJXMNPSjFb0UvDkDRHTkt9Tt4EGWe/view'),('megmod_primary.html','https://pubs.acs.org/doi/10.1021/acs.jmedchem.6c00411'),('author_guidelines.html','https://iclr.cc/Conferences/2027/AuthorGuidelines'),('ai_policy.html','https://iclr.cc/Conferences/2027/AIPolicyForAuthors')]:
 p=out/name
 if p.exists():continue
 t=time.perf_counter()
 try:
  res=requests.get(url,timeout=45);p.write_bytes(res.content);acq.append(dict(url=url,status=res.status_code,path=name,sha256=sha(p),wall_s=time.perf_counter()-t))
 except requests.RequestException as e:acq.append(dict(url=url,error=str(e),wall_s=time.perf_counter()-t))
write(out/'primary_acquisition.json',acq)
write(out/'geometry_status.json',dict(official_commit=commit,pdb_file_exists=pdb.exists(),pdb_bytes=pdb.stat().st_size if pdb.exists() else 0,rosetta_configured_executable_exists=Path(exe).exists(),RNA_FM_cached=Path('/home/shadi/.cache/torch/hub/checkpoints/RNA-FM_pretrained.pth').exists(),published_supervised_checkpoint_used=False,status='Acquisition/dependency attempt complete; inspect artifact format/coverage before declaring exact reproduction feasible.'))
