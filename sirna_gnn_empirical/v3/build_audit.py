"""Compile/render the complete paper and report actual page/float boundaries."""
from pathlib import Path
import os,sys,json,re,time,subprocess,hashlib,argparse
ROOT=Path(__file__).resolve().parents[2];P=ROOT/'papers/interaction_recoverability_iclr2027/v6';RUN=ROOT/(Path(__file__).parent/'active_run.txt').read_text().strip();D=ROOT/'docs/sirna_gnn_empirical/v3'
parser=argparse.ArgumentParser();parser.add_argument('--final',action='store_true');parser.add_argument('--no-compile',action='store_true');args=parser.parse_args()
if not args.no_compile:
 t=time.perf_counter();env=dict(os.environ,PATH='/home/shadi/.local/bin:'+os.environ['PATH']);proc=subprocess.run(['bash','build.sh'],cwd=P,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True);stamp=time.time_ns();log=RUN/'logs'/f'manuscript-build-{stamp}.log';log.write_text(proc.stdout);record=dict(command=['bash','build.sh'],cwd=str(P.relative_to(ROOT)),exit_code=proc.returncode,wall_s=time.perf_counter()-t,log=str(log.relative_to(RUN)),administrative=True)
 with (RUN/'administrative_commands.jsonl').open('a') as f:f.write(json.dumps(record)+'\n')
 print(proc.stdout[-6500:]);assert proc.returncode==0,'Manuscript compile failed; log retained.'
sys.path.insert(0,str(P/'_renderdeps'));import pymupdf
from PIL import Image,ImageOps,ImageDraw
pdf=pymupdf.open(P/'main.pdf');aux=(P/'main.aux').read_text();labels={m.group(1):int(m.group(2)) for m in re.finditer(r'\\newlabel\{([^}]+)\}\{\{[^}]*\}\{(\d+)\}',aux)};log=(P/'main.log').read_text(errors='replace');out=P/'rendered';out.mkdir(exist_ok=True);pageinfo=[];rendered=[]
for i,page in enumerate(pdf):
 pix=page.get_pixmap(matrix=pymupdf.Matrix(1.35,1.35));path=out/f'page-{i+1:03d}.png';pix.save(path);rendered.append(path);texts=page.get_text();(out/f'page-{i+1:03d}.txt').write_text(texts);blocks=[b for b in page.get_text('blocks') if b[6]==0 and len(b[4].strip())>2 and 'Under review as a conference paper' not in b[4] and not all(t.isdigit() for t in b[4].split())];pageinfo.append(dict(page=i+1,width=page.rect.width,height=page.rect.height,body_text_top=min((b[1] for b in blocks),default=None),body_text_bottom=max((b[3] for b in blocks),default=None),text_characters=len(texts)))
# Layout contact sheets show every page; main pages/figures also receive full-resolution inspection.
sheets=[]
for start in range(0,len(rendered),9):
 canvas=Image.new('RGB',(1200,1740),'#dedede');draw=ImageDraw.Draw(canvas)
 for j,path in enumerate(rendered[start:start+9]):
  im=Image.open(path).convert('RGB');im.thumbnail((390,548));x=(j%3)*400+(400-im.width)//2;y=(j//3)*580+25;canvas.paste(im,(x,y));draw.text(((j%3)*400+8,(j//3)*580+5),f'Page {start+j+1}',fill='black')
 name=out/f'contact-{start+1:03d}-{min(start+9,len(pdf)):03d}.png';canvas.save(name);sheets.append(str(name.relative_to(ROOT)))
statements=next((i+1 for i,p in enumerate(pdf) if 'ai use statement' in p.get_text().lower()),None);refs=labels.get('references:start');appendix=labels.get('appendix:start');report=dict(total_pages=len(pdf),main_pages=[1,statements-1] if statements else None,statement_pages=[statements,refs-1] if statements and refs else None,reference_pages=[refs,appendix-1] if refs and appendix else None,appendix_pages=[appendix,len(pdf)] if appendix else None,appendix_page_count=len(pdf)-appendix+1 if appendix else None,workflow_page=labels.get('fig:workflow'),main_float_pages={k:v for k,v in labels.items() if k in ['fig:workflow','fig:support','fig:activity','fig:pair','tab:data','tab:activity','tab:comparisons','tab:pairranking','tab:B3','tab:optimization','tab:calibration']},bibliography_entries=(P/'main.bbl').read_text().count('\\bibitem'),labels=labels,overfull_warnings=re.findall(r'Overfull \\[hv]box[^\n]*',log),undefined_warnings=[line for line in log.splitlines() if ('undefined' in line.lower() or 'Missing character' in line)],page_info=pageinfo,contact_sheets=sheets,pdf_sha256=hashlib.sha256((P/'main.pdf').read_bytes()).hexdigest(),visual_inspection='Pending explicit all-page visual review; rendering alone is not inspection.')
# Source placement and standalone-artwork checks supplement rendered inspection.
caption_checks=[]
for file in [P/'sections/main_text.tex',*sorted((P/'tables').glob('main_*.tex'))]:
 src=file.read_text()
 for match in re.finditer(r'\\begin\{(figure|table)\}.*?\\end\{\1\}',src,re.S):
  body=match.group();kind=match.group(1);cap=body.find(r'\caption');artifact=body.rfind(r'\includegraphics') if kind=='figure' else body.rfind(r'\end{tabular}')
  caption_checks.append(dict(file=str(file.relative_to(P)),kind=kind,caption_after_artifact=cap>artifact>=0))
report['main_caption_source_checks']=caption_checks;report['official_style_unchanged']=hashlib.sha256((P/'iclr2027_conference.sty').read_bytes()).hexdigest()==hashlib.sha256((ROOT/'papers/interaction_recoverability_iclr2027/v5/iclr2027_conference.sty').read_bytes()).hexdigest()
artwork=pymupdf.open(P/'figures/workflow.pdf');txt=' '.join(p.get_text() for p in artwork);report['workflow_artwork_has_no_embedded_caption']=not re.search(r'Figure\s+1|source.reconciliation',txt,re.I)
(D/'build_audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k not in ['labels','page_info','contact_sheets']},indent=2))
if args.final:
 assert report['main_pages']==[1,9],report['main_pages'];assert report['workflow_page']==2;assert not report['undefined_warnings'];assert not report['overfull_warnings'];assert all(v<=9 for v in report['main_float_pages'].values());assert report['bibliography_entries']>=40;assert len(report['main_float_pages'])==11;assert report['official_style_unchanged'];assert report['workflow_artwork_has_no_embedded_caption'];assert len(caption_checks)==11 and all(c['caption_after_artifact'] for c in caption_checks)
