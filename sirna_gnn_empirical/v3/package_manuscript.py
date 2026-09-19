"""Package and independently compile the complete editable manuscript."""
from common import *
import sys,subprocess,zipfile,resource,shutil
P=ROOT/'papers/interaction_recoverability_iclr2027/v6';D=ROOT/'docs/sirna_gnn_empirical/v3';sys.path.insert(0,str(P/'_renderdeps'));import pymupdf
build=json.loads((D/'build_audit.json').read_text());assert build['main_pages']==[1,9] and build['workflow_page']==2 and not build['overfull_warnings'] and not build['undefined_warnings'];assert sha(P/'main.pdf')==build['pdf_sha256']
files=[]
for p in P.rglob('*'):
 if not p.is_file():continue
 rel=p.relative_to(P)
 if any(x in rel.parts for x in ['_renderdeps','rendered','preserved_theory']):continue
 if p.suffix in ['.tex','.bib','.sty','.bst'] or rel==Path('build.sh') or (rel.parts[0]=='figures' and p.suffix=='.pdf') or (rel.parts[0]=='tables' and p.suffix=='.csv') or (rel.parts[0]=='figure_source' and p.suffix in ['.py','.md','.json']):files.append(p)
readme='Editable current empirical manuscript. Run bash build.sh with Tectonic 0.15.0. All TeX/style/bibliography/figure dependencies are included. The official style is unchanged. The source package contains the current unrestricted appendix, not reprints of superseded papers. Main text must be pages 1–9; workflow page 2. Scientific fitting code/checkpoints are in the separate code_results.zip. Figures have separate manuscript captions.\n'
(P/'SOURCE_README.txt').write_text(readme);files.append(P/'SOURCE_README.txt');manifest={str(p.relative_to(P)):dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(files)};archive=P/'manuscript_source.zip'
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
 for p in sorted(files):z.write(p,str(p.relative_to(P)))
 z.writestr('SOURCE_MANIFEST.json',json.dumps(manifest,indent=2)+'\n')
out=RUN/'package_verification/manuscript_source';out.mkdir(parents=True,exist_ok=True)
with zipfile.ZipFile(archive) as z:
 assert z.testzip() is None;z.extractall(out)
for rel,info in manifest.items():assert sha(out/rel)==info['sha256']
t=time.perf_counter();cpu=resource.getrusage(resource.RUSAGE_CHILDREN);env=dict(os.environ,PATH='/home/shadi/.local/bin:'+os.environ['PATH']);p=subprocess.run(['bash','build.sh'],cwd=out,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True);c=resource.getrusage(resource.RUSAGE_CHILDREN);log=RUN/'logs/independent-source-compile.log';log.write_text(p.stdout);assert p.returncode==0,p.stdout[-6000:]
a=pymupdf.open(P/'main.pdf');b=pymupdf.open(out/'main.pdf');assert len(a)==len(b);pagechecks=[]
for i in range(len(a)):
 textmatch=a[i].get_text()==b[i].get_text();pa=a[i].get_pixmap(matrix=pymupdf.Matrix(.6,.6));pb=b[i].get_pixmap(matrix=pymupdf.Matrix(.6,.6));imagematch=pa.width==pb.width and pa.height==pb.height and hashlib.sha256(pa.samples).digest()==hashlib.sha256(pb.samples).digest();assert textmatch and imagematch,i+1;pagechecks.append(dict(page=i+1,text_identical=textmatch,render_identical=imagematch))
report=dict(archive=str(archive.relative_to(ROOT)),archive_sha256=sha(archive),source_files=len(files),independent_command=['bash','build.sh'],independent_cwd=str(out.relative_to(ROOT)),exit_code=p.returncode,wall_s=time.perf_counter()-t,child_cpu_s=c.ru_utime+c.ru_stime-cpu.ru_utime-cpu.ru_stime,log=str(log.relative_to(RUN)),pages=len(a),all_page_text_and_renders_identical=True,page_checks=pagechecks,main_pages=build['main_pages'],reference_pages=build['reference_pages'],appendix_pages=build['appendix_pages'])
write(RUN/'package_verification/manuscript_source.json',report);write(D/'source_zip_verification.json',report)
# Caption-free artwork is separate from its descriptive LaTeX caption and scientific specification.
wfiles=[P/'figures'/f'workflow.{ext}' for ext in ['pdf','svg','png']]+[P/'figures/workflow_caption.tex']+[p for p in (P/'figure_source').iterdir() if p.is_file() and p.suffix in ['.py','.md','.json']]
with zipfile.ZipFile(P/'workflow_handoff.zip','w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in wfiles:z.write(p,str(p.relative_to(P)))
 z.writestr('README.txt','The PDF/SVG/PNG contain no embedded Figure 1 caption. The separate LaTeX caption is placed below the artwork in the paper. WORKFLOW_SPEC.md defines the actual forward pass, controls and separate activity/pair objectives. The editable vector generator uses Matplotlib/CairoSVG as recorded; no model fitting is needed.\n')
 z.writestr('MANIFEST.json',json.dumps({str(p.relative_to(P)):sha(p) for p in wfiles},indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='page_checks'},indent=2))
