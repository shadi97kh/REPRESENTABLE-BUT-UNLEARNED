"""Administrative compile/check of the edited manuscript; never regenerate accepted prose."""
from core import *
import subprocess,re,shutil,zipfile
P=ROOT/'papers/interaction_recoverability_iclr2027/v3';D=ROOT/'docs/sirna_gnn_empirical/v2';D.mkdir(exist_ok=True)
result=subprocess.run(['/home/shadi/.local/bin/tectonic','--keep-logs','--keep-intermediates','main.tex'],cwd=P,text=True,capture_output=True);(P/'private/final_build.log').write_text(result.stdout+result.stderr);assert result.returncode==0
aux=(P/'main.aux').read_text();label=lambda k:int(re.search(r'\\newlabel\{'+re.escape(k)+r'\}\{\{[^}]*\}\{(\d+)\}',aux)[1]);assert label('main:end')==9
log=(P/'main.log').read_text();assert not re.search(r'(?:Reference|Citation).*undefined|multiply defined|Overfull \\[hv]box',log), 'unresolved citation/reference/overflow'
info=subprocess.check_output(['pdfinfo',str(P/'main.pdf')],text=True);pages=int(re.search(r'Pages:\s+(\d+)',info)[1]);audit=dict(main_pages=[1,9],statements=[10,10],references=[label('references:start'),label('references:end')],appendix=[label('appendix:start'),pages],appendix_pages=pages-label('appendix:start')+1,total_pages=pages,citations=len(re.findall(r'\\bibitem', (P/'main.bbl').read_text())),official_style_sha256=sha(P/'iclr2027_conference.sty'),preserved_style_sha256=sha(P.parent/'v2/iclr2027_conference.sty'),undefined_references=False,overflow=False,pdf_sha256=sha(P/'main.pdf'))
assert audit['official_style_sha256']==audit['preserved_style_sha256'];write(RUN/'paper_build.json',audit);write(D/'build_audit.json',audit)
# A results-focused extract includes current main text and the complete current follow-up appendix; historical proofs remain in the full PDF.
start=label('appendix:start');hist=label('app:historical');temp=P/'private/extract';temp.mkdir(exist_ok=True)
subprocess.run(['pdfseparate',str(P/'main.pdf'),str(temp/'page-%d.pdf')],check=True)
indices=list(range(1,10))+list(range(start,hist));subprocess.run(['pdfunite',*[str(temp/f'page-{i}.pdf') for i in indices],str(P/'results_extract.pdf')],check=True)
write(D/'results_extract_mapping.json',dict(full_pdf_pages=indices,excluded_statements_and_references=True,historical_complete_evidence_retained_in_full_pdf=True))
print(audit)
