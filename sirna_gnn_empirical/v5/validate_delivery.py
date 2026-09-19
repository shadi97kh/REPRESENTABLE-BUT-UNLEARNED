"""External build, exact source/ZIP validation, rendered-page inventory and clean ZIP build."""
from common import *
import resource,shlex,zipfile,tempfile,re,shutil,sys
HERE=Path(__file__).resolve().parent
import subprocess,resource,shlex,zipfile,tempfile,re,shutil,sys
sys.path.insert(0,str(ROOT/'papers/interaction_recoverability_iclr2027/v6/_renderdeps'))
import pymupdf as fitz
from PIL import Image,ImageOps,ImageDraw
P=ROOT/'papers/interaction_recoverability_iclr2027/v10';S=P/'source';A=P/'artifacts';D=ROOT/'docs/sirna_gnn_empirical/v5';A.mkdir(exist_ok=True);D.mkdir(exist_ok=True);logs=A/'validation_logs';logs.mkdir(exist_ok=True);records=[]
def command(cmd,name,cwd=ROOT):
 t=time.monotonic();c=resource.getrusage(resource.RUSAGE_CHILDREN);result=subprocess.run(cmd,cwd=cwd,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);z=resource.getrusage(resource.RUSAGE_CHILDREN);(logs/(name+'.log')).write_text(result.stdout);records.append(dict(command=shlex.join(map(str,cmd)),exit_code=result.returncode,wall_s=time.monotonic()-t,child_cpu_s=z.ru_utime+z.ru_stime-c.ru_utime-c.ru_stime,max_rss_kib=z.ru_maxrss,log=str((logs/(name+'.log')).relative_to(ROOT))));write(A/'validation_commands.json',records)
 if result.returncode:print(result.stdout[-5000:]);raise RuntimeError(name)
command([sys.executable,str(HERE/'validate_minimal_tex.py'),str(S),'--report',str(A/'minimal_source_validation.json')],'minimal_source')
command(['bash',str(P/'typesetting/build.sh'),str(S),str(A/'build')],'external_build')
log=(A/'build/main.log').read_text();bad=re.findall(r'(?:Overfull \\[hv]box[^\n]*|[^\n]*(?:undefined references|undefined citations|multiply defined|Reference .* undefined|Citation .* undefined)[^\n]*)',log,re.I);write(A/'typesetting_warnings.json',dict(blocking_warnings=bad));assert not bad,bad[:12]
aux=(A/'build/main.aux').read_text();labels={k:int(p) for k,p in re.findall(r'\\newlabel\{([^}]+)\}\{\{[^}]*\}\{(\d+)\}',aux)};doc=fitz.open(A/'build/main.pdf');boundary=dict(main_start=1,main_end=labels['main:end'],statements_start=labels['main:end']+1,references_start=labels['references:start'],references_end=labels['references:end'],appendix_start=labels['appendix:start'],appendix_end=len(doc),appendix_pages=len(doc)-labels['appendix:start']+1,total_pages=len(doc),workflow_page=labels['fig:workflow']);write(A/'page_boundaries.json',boundary);assert boundary['main_end']==9,boundary;assert boundary['workflow_page']==2,boundary;assert len(doc[1].get_drawings())>0; assert len(doc[1].get_images())==0
shutil.copy2(A/'build/main.pdf',A/'main.pdf');first=fitz.open();first.insert_pdf(doc,from_page=0,to_page=8);first.save(A/'main_text_9_pages.pdf')
# Fully render each page; contact sheets aid whole-document visual inspection.
render=A/'rendered';render.mkdir(exist_ok=True);pageinfo=[];thumbs=[]
for i,page in enumerate(doc):
 pix=page.get_pixmap(matrix=fitz.Matrix(110/72,110/72),alpha=False);pix.save(render/f'page-{i+1:03}.png');im=Image.open(render/f'page-{i+1:03}.png').convert('RGB');im.thumbnail((540,720));tile=Image.new('RGB',(570,760),'#dddddd');tile.paste(im,((570-im.width)//2,25));ImageDraw.Draw(tile).text((12,5),f'PDF page {i+1}',fill='black');thumbs.append(tile)
 blocks=[b for b in page.get_text('blocks') if b[4].strip()];pageinfo.append(dict(page=i+1,text_characters=len(page.get_text()),images=len(page.get_images()),minimum_text_x=min([b[0] for b in blocks],default=0),maximum_text_x=max([b[2] for b in blocks],default=0)))
 if len(thumbs)==6 or i==len(doc)-1:
  sheet=Image.new('RGB',(1710,1520),'#bbb')
  for j,tile in enumerate(thumbs):sheet.paste(tile,((j%3)*570,(j//3)*760))
  sheet.save(render/f'contact-{(i//6)+1:02}.jpg',quality=88);thumbs=[]
write(A/'render_inventory.json',pageinfo)
# ZIP has exactly the source's four entries, including the one assets directory.
zipname=A/'manuscript_source.zip'
with zipfile.ZipFile(zipname,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 z.writestr('figures/','')
 for p in sorted(S.rglob('*')):
  if p.is_file():z.write(p,str(p.relative_to(S)))
command([sys.executable,str(HERE/'validate_minimal_tex.py'),str(zipname),'--report',str(A/'minimal_zip_validation.json')],'minimal_zip')
with tempfile.TemporaryDirectory(prefix='sirna-v10-independent-',dir='/tmp') as tmp:
 tmp=Path(tmp);clean=tmp/'source';clean.mkdir();deps=tmp/'typesetting';shutil.copytree(P/'typesetting',deps)
 with zipfile.ZipFile(zipname) as z:z.extractall(clean)
 command(['bash',str(deps/'build.sh'),str(clean),str(tmp/'output')],'independent_zip_build')
 independent=fitz.open(tmp/'output/main.pdf');assert len(independent)==len(doc)
 mismatches=[i+1 for i in range(len(doc)) if independent[i].get_text()!=doc[i].get_text()];assert not mismatches,mismatches
 ilog=(tmp/'output/main.log').read_text();assert not re.search(r'Overfull \\[hv]box|There were undefined references',ilog)
 shutil.copy2(tmp/'output/main.pdf',A/'independent_zip_build.pdf')
 write(A/'independent_build.json',dict(passed=True,pages=len(independent),all_page_text_equal=True,source_zip_sha256=sha(zipname),external_dependencies={p.name:sha(p) for p in deps.glob('*') if p.is_file()}))
with zipfile.ZipFile(A/'typesetting_dependencies.zip','w',zipfile.ZIP_DEFLATED) as z:
 for p in sorted((P/'typesetting').glob('*')):
  if p.is_file():z.write(p,p.name)
proofs=re.findall(r'\\begin\{proof\}.*?\\end\{proof\}',(ROOT/'papers/interaction_recoverability_iclr2027/v9/source/appendix.tex').read_text(),re.S);app=(S/'appendix.tex').read_text();assert all(p in app for p in proofs),'A retained mathematical proof changed or disappeared'
write(A/'proof_preservation.json',dict(original_complete_proofs=len(proofs),all_original_proofs_preserved_verbatim=True,current_complete_proofs=len(re.findall(r'\\begin\{proof\}',app))))
write(D/'build_audit.json',dict(boundaries=boundary,zero_undefined_references=True,zero_overfull_boxes=True,independent_zip_build=True,all_pages_rendered=True,visual_inspection='Rendered contact sheets and full main pages require human-facing agent inspection; recorded separately after inspection.',source_entries=sorted(p.name for p in S.iterdir()),official_style_sha256=sha(P/'typesetting/iclr2027_conference.sty')))
print(json.dumps(boundary,indent=2))
