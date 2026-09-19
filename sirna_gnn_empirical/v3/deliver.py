"""Final preservation, scientific/source packages and honest delivery accounting."""
from common import *
import subprocess,sys,zipfile,resource,shlex
D=ROOT/'docs/sirna_gnn_empirical/v3';P=ROOT/'papers/interaction_recoverability_iclr2027/v6';out=RUN/'paper_artifacts';out.mkdir(exist_ok=True)
build=json.loads((D/'build_audit.json').read_text());assert build['main_pages']==[1,9] and build['workflow_page']==2 and not build['overfull_warnings'] and not build['undefined_warnings'];assert build.get('visual_inspection_completed') is True,'Explicit complete-page visual inspection must precede delivery';assert sha(P/'main.pdf')==build['pdf_sha256']
t0=time.perf_counter();c0=time.process_time();operations=[];attempt_stamp=str(time.time_ns());attempt_dir=RUN/'delivery_attempts'/attempt_stamp;attempt_dir.mkdir(parents=True,exist_ok=True)
def command(script):
 t=time.perf_counter();cpu=resource.getrusage(resource.RUSAGE_CHILDREN);cmd=[sys.executable,str(HERE/script)];p=subprocess.run(cmd,cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True);c=resource.getrusage(resource.RUSAGE_CHILDREN);log=attempt_dir/(Path(script).stem+'.log');log.write_text(p.stdout);operations.append(dict(command=shlex.join(cmd),exit_code=p.returncode,wall_s=time.perf_counter()-t,child_cpu_s=c.ru_utime+c.ru_stime-cpu.ru_utime-cpu.ru_stime,log=str(log.relative_to(RUN))));write(attempt_dir/'operations.json',operations);print(p.stdout[-4000:],flush=True);assert p.returncode==0,script
command('main_cell_provenance.py');command('pair_gradient_check.py');command('audit_results.py');command('package_manuscript.py');command('write_delivery_docs.py')
base=json.loads((RUN/'preservation_before.json').read_text());changed=[];missing=[];count=0;totalbytes=0
for rel,info in base.items():
 p=ROOT/rel
 if not p.exists():missing.append(rel);continue
 if sha(p)!=info['sha256']:changed.append(rel)
 count+=1;totalbytes+=p.stat().st_size
pres=dict(inventoried_files=len(base),checked_files=count,bytes_hashed=totalbytes,changed=changed,missing=missing,all_preserved=not(changed or missing),method='SHA256 content comparison against inventory captured before current campaign; existing untracked user files included.')
write(RUN/'preservation_verification.json',pres);write(D/'integrity_audit.json',pres);assert pres['all_preserved'],pres
status=subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True);(RUN/'git_status_after.txt').write_text(status)
# All current scientific fit objects and records are included. Administrative citation downloads and unused external representation assets are omitted deliberately.
files=set();omitted=[]
for p in RUN.rglob('*'):
 if not p.is_file():continue
 rel=p.relative_to(RUN)
 if rel.parts[0] in ['paper_artifacts','package_verification'] or p.name in ['campaign.lock','deliver.done.json'] or '__pycache__' in rel.parts:continue
 if rel.parts[:2]==('sources','citations') and p.suffix not in ['.json']:
  omitted.append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p),reason='Administrative downloaded primary text; citation links/acquisition script and audit supplied, not needed for reported scientific prediction.'));continue
 if p.name in ['unimol_1b_emb_dict.pkl','cofold_results.pkl']:
  omitted.append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p),reason='Unused external representation asset from blocked exact-method attempt; public acquisition provenance supplied.'));continue
 files.add(p)
for root in [HERE,D]:
 for p in root.rglob('*'):
  if p.is_file() and not any(x in p.relative_to(root).parts for x in ['_deps','__pycache__']) and p.name!='WORK_IN_PROGRESS.md':files.add(p)
# Exact historical input dependencies for core training, legacy recomputation and source reconciliation.
oldnames=[f'{s}_{t}' for s in ['train','validation','test_internal','test_APP'] for t in ['observations.jsonl','graphs.npz']]+['eligible_pairs.jsonl','eligible_rectangles.jsonl','graph_manifest.json','heldout_predictions.csv','chemical_effect_metrics.csv']
prevnames=['protocol.json','grouped_membership.csv','b2_membership.csv','all_B1_predictions.csv','grouped_predictions.csv','external_predictions.csv','external_observations.jsonl','external_graphs.npz','b2_predictions.csv','relation_support.csv','data_gradients.csv','frozen_replay.json','source_audit/davis_original_supplement.xlsx','source_audit/davis_orientation_rows.csv']
for parent,names in [(OLD,oldnames),(PREV,prevnames)]:
 for name in names:
  p=parent/name;assert p.exists(),p;files.add(p)
# Source-data attribution/acquisition and graph builders remain inspectable without old fitted checkpoints.
for version in ['v1','v2']:
 for p in (ROOT/'sirna_gnn_empirical'/version).rglob('*'):
  if p.is_file() and p.suffix in ['.py','.sh','.json','.md','.txt'] and not any(x in p.parts for x in ['__pycache__','_deps']):files.add(p)
for name in ['Codex_Paper_Revision_and_Experiments.md','Reviewer_Critique.md','Review_Integration_Notes.md']:
 p=ROOT/'handoffs/codex_revision_v1'/name
 if p.exists():files.add(p)
# Keep official code acquisition reproducible at pinned commits; no fitted efficacy checkpoint is substituted.
readme='''# Current scientific code/results package

This archive uses repository-relative paths. It contains every current fit checkpoint, optimizer state, history, prediction, frozen protocol, metric and plotted mark, together with the exact admitted observation/graph inputs and historical prediction inputs used by current analysis. The original complete local run remains separately preserved.

The current manuscript is a separate manuscript_source.zip. Extract it into papers/interaction_recoverability_iclr2027/v6/ before manuscript stages. Install the environment recorded in the run's environment.json. The original executed entry point is bash sirna_gnn_empirical/v3/run_all.sh; its interpreter is /opt/miniforge3/bin/python. No scientific fit needs to be rerun to inspect the reported evidence.

Use VERIFY_PACKAGE.py here to verify every included file. Every original per-fit hash is preserved and all referenced fit files are included. Original .done records also cover administrative downloads produced concurrently during training. This focused archive intentionally omits downloaded citation full texts and unused published-representation assets; therefore full original-stage .done verification requires the complete local archive, not only this focused package. PACKAGE_OMISSIONS.json identifies every omitted current-run administrative source and its original hash. This is not a claim that omitted files were absent during execution.

For a fresh scientific rerun, work in a copy, select a NEW run directory through sirna_gnn_empirical/v3/active_run.txt, copy preservation_before.json and the admitted primary B3 protocol/data and primary-source records into that directory, and use the recorded stage commands. Do not delete the completed original run or overwrite historical inputs. The fresh protocol uses the packaged prior grouped membership and exact observation/graph files. The source stage additionally expects author code at the pinned paths below; clone those official commits before rerunning source checks. That stage records the known reproduction blockers, not a fitted published predictor.

Official code:
  git clone https://github.com/tanwenchong/ENsiRNA handoffs/codex_revision_v1/ENsiRNA_official
  git -C handoffs/codex_revision_v1/ENsiRNA_official checkout 028824341635903f3c661f5d1cc737de106493d5
  git clone https://github.com/YuantingChen111/MEG-mod.git handoffs/codex_revision_v1/MEG_mod_official
  git -C handoffs/codex_revision_v1/MEG_mod_official checkout c335a4c69d56ef73754677da8bfe627c8112b350
These acquisition commands are reproduction instructions; their original clone/inspection attempts are recorded separately. Exact ENsiRNA geometry and the MEG modification-token path remain blockers. No approximation replaces them.

Historical code files provide parser/graph provenance; historical runs and checkpoints are not presented as new fits. Public primary measurement source links, original hashes and normalization details remain in the current source records and manuscript. APP/S7 evaluation is retrospective. The separate primary B3 panel has one background/shared reference. No new wet-lab or causal validation is claimed.
'''
(out/'README_code_results.md').write_text(readme);write(out/'PACKAGE_OMISSIONS.json',omitted)
verify='''from pathlib import Path
import json,hashlib
root=Path(__file__).resolve().parent
m=json.loads((root/'PACKAGE_MANIFEST.json').read_text())
for rel,info in m.items():
 p=root/rel
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 assert h.hexdigest()==info['sha256'],rel
print('Verified',len(m),'included files')
'''
manifest={str(p.relative_to(ROOT)):dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(files)};archive=out/'code_results.zip'
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=3,allowZip64=True) as z:
 for p in sorted(files):z.write(p,str(p.relative_to(ROOT)))
 z.writestr('README.md',readme);z.writestr('VERIFY_PACKAGE.py',verify);z.writestr('PACKAGE_MANIFEST.json',json.dumps(manifest,indent=2)+'\n');z.writestr('PACKAGE_OMISSIONS.json',json.dumps(omitted,indent=2)+'\n')
# Stream and hash every archived file without duplicating multi-GB checkpoints on disk.
with zipfile.ZipFile(archive) as z:
 assert z.testzip() is None
 for rel,info in manifest.items():
  h=hashlib.sha256()
  with z.open(rel) as f:
   for data in iter(lambda:f.read(4*1024*1024),b''):h.update(data)
  assert h.hexdigest()==info['sha256'],rel
report=dict(archive=str(archive.relative_to(ROOT)),sha256=sha(archive),bytes=archive.stat().st_size,included_files=len(manifest),uncompressed_bytes=sum(i['bytes'] for i in manifest.values()),every_included_file_stream_hash_verified=True,all_current_fit_artifacts_included=True,omitted_administrative_files=len(omitted),source_zip_separately_verified=True,preservation=pres,operations=operations,administrative_process_cpu_s=time.process_time()-c0,administrative_wall_s=time.perf_counter()-t0,qualification='Packaging costs are administrative and overlap runner deliver-stage time. Interactive inspection/authoring/browsing and prior parent hashing retain additional unmetered nonzero costs. ZIP snapshot predates its own completion record to avoid recursive self-hashes.')
write(out/'package_verification.json',report);write(D/'package_verification.json',report)
review=out/'review_packet.zip'
reviewfiles=[P/'main.pdf',P/'manuscript_source.zip',P/'workflow_handoff.zip',RUN/'protocol.md',RUN/'protocol.json',RUN/'experiment_matrix.csv',RUN/'resource_accounting.json',RUN/'main_cell_provenance.json',RUN/'figure_provenance.json',RUN/'final_scientific_checks.json',RUN/'preservation_verification.json',out/'package_verification.json']+list(D.glob('*.md'))+list(D.glob('*.json'))+list((RUN/'tables').glob('*.csv'))
reviewfiles=[p for p in reviewfiles if p.name!='WORK_IN_PROGRESS.md']
request='Draft review request; not sent. Please assess whether the controlled support/message comparison, training-fitted constants, conditional ensemble uncertainty, separate ten-seed variability, B2 mean-control and within-assay diagnostics, and single-background primary B3 evaluation support the stated claims. Check the primary/released-label discrepancies, blocked exact published methods and retrospective APP/S7 scope. No GNN novelty, causal mechanism, private-anchor theorem coverage or new wet-lab validation is claimed. The compact packet contains the full paper/sources/audits/tables. All checkpoint/optimizer/prediction artifacts are in the separately supplied code_results.zip; the compact packet does not replace that scientific archive. No favorable result is assumed.\n'
with zipfile.ZipFile(review,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
 for p in sorted(set(reviewfiles)):z.write(p,str(p.relative_to(ROOT)))
 z.writestr('REVIEW_REQUEST.md',request);z.writestr('REVIEW_MANIFEST.json',json.dumps({str(p.relative_to(ROOT)):sha(p) for p in sorted(set(reviewfiles))},indent=2)+'\n')
with zipfile.ZipFile(review) as z:assert z.testzip() is None
artifacts={str(p.relative_to(ROOT)):dict(sha256=sha(p),bytes=p.stat().st_size) for p in [P/'main.pdf',P/'manuscript_source.zip',P/'workflow_handoff.zip',P/'figures/workflow.pdf',P/'figures/workflow.svg',P/'figures/workflow.png',archive,review]};write(RUN/'delivery_manifest.json',dict(artifacts=artifacts,build=build,preservation=pres,package_verification=str((out/'package_verification.json').relative_to(ROOT)),complete_core_campaign=True,published_exact_reproduction='BLOCKED, zero fits',wetlab_validation=False,push_performed=False));write(D/'delivery_manifest.json',json.loads((RUN/'delivery_manifest.json').read_text()))
print('Delivered verified packages:',json.dumps(artifacts,indent=2),flush=True)
