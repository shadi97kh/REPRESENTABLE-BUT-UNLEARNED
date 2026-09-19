"""Compile the complete manuscript and package existing scientific evidence. No new fits."""
import argparse,datetime,hashlib,json,re,shutil,subprocess,sys,time,zipfile
from pathlib import Path
from common import write_json
ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);r=ap.parse_args().run
code=Path(__file__).resolve().parent;repo=code.parents[1];base=repo/'papers/interaction_recoverability_iclr2027/v2';active=(repo/(code/'active_run.txt').read_text().strip()).resolve();p=base if r.resolve()==active else r/'manuscript'
assert (r/'figures.done.json').exists()
if p!=base and not p.exists():
 shutil.copytree(base,p,ignore=shutil.ignore_patterns('private','*.zip','*.aux','*.log','*.out','*.bbl','*.blg','main.pdf'))
(p/'private').mkdir(exist_ok=True)
for f in (r/'figures').glob('*.pdf'):shutil.copy2(f,p/'figures'/f.name)
subprocess.run([sys.executable,str(base/'scripts/generate_support.py'),'--paper-dir',str(p),'--repo',str(repo),'--run',str(r)],check=True)
t0=time.monotonic();subprocess.run(['bash','build.sh'],cwd=p,check=True)
subprocess.run(['pdftotext','-layout','main.pdf','private/rendered_text.txt'],cwd=p,check=True)
aux=(p/'main.aux').read_text();log=(p/'main.log').read_text();txt=(p/'private/rendered_text.txt').read_text();pages=txt.split('\f');pages=pages[:-1] if not pages[-1].strip() else pages
labels={m[1]:{'number':m[2],'page':int(m[3])} for m in re.finditer(r'\\newlabel\{([^}]+)\}\{\{([^}]*)\}\{(\d+)\}',aux)}
assert labels['main:end']['page']==9,labels['main:end'];assert labels['references:start']['page']>9
bad=[x for x in log.splitlines() if any(z in x for z in ['undefined references','Citation `','multiply defined','Overfull \\hbox','Overfull \\vbox','Missing character:'])];assert not bad,bad
src='\n'.join(f.read_text() for f in p.rglob('*.tex') if 'private' not in f.parts and 'workflow_handoff' not in f.parts);refs=re.findall(r'\\(?:eqref|ref|pageref)\{([^}]+)\}',src);assert all(x in labels for x in refs),sorted(set(refs)-set(labels))
bib=(p/'main.bbl').read_text();bibcount=bib.count('\\bibitem');assert bibcount>=42
main=(p/'sections/main_text.tex').read_text();assert main.count('\\begin{figure}')>=4;assert not (p/'figures/workflow.pdf').exists(),'This stage must preserve the undrawn workflow specification'
assert '/home/' not in txt and 'shadi' not in txt.lower()
style=json.loads((base/'style_origin.json').read_text());assert all(hashlib.sha256((p/n).read_bytes()).hexdigest()==h for n,h in style['files'].items())
report={'compiled_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'main_pages':[1,9],'excluded_statements_pages':[10,labels['references:start']['page']-1],'reference_pages':[labels['references:start']['page'],labels['references:end']['page']],'appendix_pages':[labels['app:roadmap']['page'],labels['appendix:end']['page']],'appendix_page_count':labels['appendix:end']['page']-labels['app:roadmap']['page']+1,'total_pages':len(pages),'bibliography_entries':bibcount,'research_method_software_papers':bibcount-2,'assay_corrections':2,'main_figure_slots':main.count('\\begin{figure}'),'main_tables':main.count('\\begin{table}'),'appendix_figures':(p/'appendices/results.tex').read_text().count('\\begin{figure}'),'undefined_references':0,'overfull_boxes':0,'official_style_unchanged':True,'workflow_drawn':False,'main_page_word_counts':[len(x.split()) for x in pages[:9]],'all_labels':labels,'build_and_checks_wall_s':time.monotonic()-t0,'visual_inspection':'Main pages and appendix renders reviewed separately in visual_inspection.md; no assertion of external expert review.'}
write_json(p/'private/final_build_audit.json',report)
art=r/'paper_artifacts';art.mkdir(exist_ok=True)
shutil.copy2(p/'main.pdf',art/'main.pdf');write_json(art/'final_build_audit.json',report)
for doc in ['citation_audit.md','claim_audit.md','appendix_audit.md']:
 shutil.copy2(repo/'docs/sirna_gnn_empirical/v1'/doc,art/doc)
with zipfile.ZipFile(art/'workflow_handoff.zip','w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for f in sorted((p/'workflow_handoff').iterdir()):
  if f.is_file():z.write(f,'workflow_handoff/'+f.name)
with zipfile.ZipFile(art/'manuscript_source.zip','w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for f in sorted(p.rglob('*')):
  rel=f.relative_to(p)
  if not f.is_file() or any(x in rel.parts for x in ['private','__pycache__']):continue
  if f.suffix not in ['.tex','.bib','.sty','.bst','.sh','.py','.md','.json','.pdf']:continue
  if f.name=='main.pdf':continue
  z.write(f,'manuscript/'+str(rel))
# Anonymous complementary scientific supplement: exact scientific files, without machine/user logs.
excluded={'hardware_inventory.json','preservation_before.json','git_status_before.txt','commands.jsonl','campaign.lock'}
with zipfile.ZipFile(art/'code_results.zip','w',zipfile.ZIP_DEFLATED,compresslevel=1) as z:
 for f in sorted(code.iterdir()):
  if f.is_file() and f.suffix in ['.py','.sh','.txt','.md','.json']:z.write(f,'sirna_gnn_empirical/v1/'+f.name)
 for name in ['siRNA_Transfer_Source_Evidence.json'] :z.write(repo/name,name)
 for f in sorted(p.rglob('*')):
  rel=f.relative_to(p)
  if f.is_file() and not any(x in rel.parts for x in ['private','__pycache__']) and f.suffix in ['.tex','.bib','.sty','.bst','.sh','.py','.md','.json','.pdf'] and f.name!='main.pdf':z.write(f,'papers/interaction_recoverability_iclr2027/v2/'+str(rel))
 for f in sorted((repo/'docs/sirna_gnn_empirical/v1').glob('*.md')):z.write(f,'docs/sirna_gnn_empirical/v1/'+f.name)
 manifest=[]
 for f in sorted(r.rglob('*')):
  rel=f.relative_to(r)
  if not f.is_file() or any(x in rel.parts for x in ['paper_artifacts','manuscript','logs','superseded']):continue
  if f.name in excluded or any(t in f.name for t in ['.done.json','.running.json','.failed-']):continue
  raw=f.read_bytes()
  if b'/home/shadi' in raw:continue
  target='runs/sirna_gnn_empirical/'+r.name+'/'+str(rel);z.writestr(target,raw);manifest.append({'path':target,'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw)})
 z.writestr('scientific_artifact_manifest.json',json.dumps(manifest,indent=2)+'\n');z.writestr('SUPPLEMENT_README.md','Scientific code, acquired raw inputs, reconciled observations, frozen splits, features, fit checkpoints/histories, committed predictions, results and figures. Machine-specific original logs and preservation records remain in the local evidence run and are omitted here for anonymity. Stage completion markers are omitted: use a NEW --run path for a full replay; do not overwrite the supplied frozen scientific evidence. The complete procedure and executed commands are in the manuscript appendix. No new wet-lab validation is claimed.\n')
for name in ['manuscript_source.zip','workflow_handoff.zip']:shutil.copy2(art/name,p/name)
write_json(art/'artifact_hashes.json',{f.name:{'sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'bytes':f.stat().st_size} for f in art.iterdir() if f.is_file() and f.name!='artifact_hashes.json'})
write_json(r/'paper.outputs.json',[str(f.relative_to(r)) for f in sorted(art.iterdir()) if f.is_file()]);print(json.dumps(report,indent=2))
