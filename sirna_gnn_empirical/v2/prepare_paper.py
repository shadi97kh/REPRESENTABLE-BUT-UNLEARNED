from core import *
import re,shutil
P=ROOT/'papers/interaction_recoverability_iclr2027/v3';V2=P.parent/'v2';V1=P.parent/'v1'
for pat in ['*.sty','*.bst','*.bib','style_origin.json','build.sh']:
 for p in V2.glob(pat):shutil.copy2(p,P/p.name)
for p in (RUN/'figures').glob('*.pdf'):shutil.copy2(p,P/'figures'/p.name)
# Merge historical bibliography keys without replacing accepted records.
existing=''.join((P/f).read_text() for f in ['references.bib','references_new.bib']);keys=set(re.findall(r'@\w+\s*\{\s*([^,]+),',existing));bib=(V1/'references.bib').read_text();entries=re.split(r'(?=@\w+\s*\{)',bib);new=[]
for entry in entries:
 m=re.match(r'@\w+\s*\{\s*([^,]+),',entry)
 if m and m[1] not in keys:new.append(entry);keys.add(m[1])
(P/'references_theory.bib').write_text(''.join(new))
def bottom_captions(s):
 # Balanced brace scanning preserves multiline captions and labels exactly.
 for env in ['longtable','table','figure']:
  pattern=re.compile(r'\\begin\{'+env+r'\}.*?\\end\{'+env+r'\}',re.S)
  def fix(m):
   text=m[0];start=text.find('\\caption{')
   if start<0:return text
   j=start+len('\\caption');depth=0;end=j
   for pos in range(j,len(text)):
    if text[pos]=='{' and (pos==0 or text[pos-1]!='\\'):depth+=1
    if text[pos]=='}' and text[pos-1]!='\\':
     depth-=1
     if depth==0:end=pos+1;break
   match=re.match(r'\s*\\label\{[^}]+\}',text[end:])
   if match:end+=match.end()
   cap=text[start:end];tail=text[end:]
   if env=='longtable' and tail.startswith('\\\\'):tail=tail[2:]
   body=text[:start]+tail;ep=body.rfind('\\end{'+env+'}');return body[:ep]+'\n'+cap+('\\\\\n' if env=='longtable' else '\n')+body[ep:]
  s=pattern.sub(fix,s)
 return s
for tag,src in [('historical',V2),('theory',V1)]:
 dst=P/tag;dst.mkdir(exist_ok=True)
 for sub in ['sections','appendices','figures','workflow_handoff']:
  if not (src/sub).exists():continue
  for p in (src/sub).rglob('*'):
   if not p.is_file() or p.suffix not in ['.tex','.pdf','.png','.jpg','.json','.md']:continue
   target=dst/p.relative_to(src);target.parent.mkdir(parents=True,exist_ok=True)
   if p.suffix=='.tex':
    s=p.read_text();s=re.sub(r'\\(label|ref|eqref|pageref)\{([^}]+)\}',lambda m:'\\'+m[1]+'{'+tag+':'+m[2]+'}',s)
    s=re.sub(r'\\(input|includegraphics|IfFileExists)(\[[^\]]*\])?\{([^}]+)\}',lambda m:'\\'+m[1]+(m[2] or '')+'{'+tag+'/'+m[3]+'}',s)
    # Historical main-text pagination is unnecessary in the unlimited appendix.
    s=s.replace('\\clearpage','').replace('\\begin{abstract}','\\paragraph{Historical abstract.}').replace('\\end{abstract}','')
    target.write_text(bottom_captions(s))
   else:shutil.copy2(p,target)
# The final workflow is reserved once in the current main text; preserve historical scientific content without duplicate empty artwork.
for tag in ['historical','theory']:
 p=P/tag/'sections/main_text.tex';s=p.read_text();s=re.sub(r'\\begin\{figure\}.*?\\end\{figure\}',lambda m:('Historical workflow specification is superseded by the current reserved workflow slot.\\label{'+tag+':fig:workflow}\n') if 'workflow' in m[0] else m[0],s,flags=re.S);p.write_text(s)
main=V1.joinpath('main.tex').read_text().split('\\begin{document}')[0]
main=main.replace('Interaction Extrapolation from Five Source Laws: Attainment and Exponential Conditioning','Support-Aware siRNA Prediction: A Targeted Measured-Data Follow-up')
main=re.sub(r'\\title\{.*?\}\n\\author',r'\\title{Support-Aware siRNA Prediction:\\\\Measured Activity, Grouped Validation,\\\\and Supervised Chemistry Effects}\n\\author',main,flags=re.S)
main+=r'''
\newcommand{\MSE}{\operatorname{MSE}}
\begin{document}
\maketitle
\input{sections/main_text}
\label{main:end}\clearpage
\input{sections/statements}\clearpage
\label{references:start}
\bibliography{references,references_new,references_theory}
\bibliographystyle{iclr2027_conference}
\label{references:end}\clearpage
\appendix\numberwithin{equation}{section}
\label{appendix:start}
\input{appendices/followup}
\clearpage
\section{Preserved original empirical campaign: historical record}\label{app:historical}
The following protocol, results and limitations describe the completed v1 empirical campaign, before the targeted follow-up. Statements about unimplemented follow-up work or unadmitted Davis data are historical. The current results and source-specific admission decisions appear above. Original numerical outcomes are retained, not replaced. Figure/table captions have been moved below their content as accepted in the review; no historical numerical record is changed.
\input{historical/sections/main_text}
'''
for name in ['roadmap','data','model','training','evaluation','results','reproduction','scope']:main+='\\input{historical/appendices/'+name+'}\n'
main+=r'''
\clearpage
\section{Preserved calibrated five-law theory and complete proofs}\label{app:theory}
This independent mathematical experiment is reproduced to make the theoretical claims self-contained. Its historical pilot and reporter comparisons are distinct from both empirical campaigns above. The theorem does not cover either learned GNN. Original assumptions, constructions, proofs and negative qualifications remain in force. Historical workflow artwork is omitted because its final scientific specification belongs to the current handoff.
\input{theory/sections/main_text}
'''
for name in ['roadmap','notation','geometry','estimation','conditioning','hard_profiles','hard_score','lower_bound','expanded_arguments','corollary','prior','evidence','application','visual_evidence']:main+='\\input{theory/appendices/'+name+'}\n'
main+='\\label{appendix:end}\n\\end{document}\n';(P/'main.tex').write_text(main)
(P/'sections/statements.tex').write_text(r'''\section*{Reproducibility statement}
The complete appendix specifies the measured-data inputs, source checks, grouped partitions, model equations, training schedules, full fit grid and every reported endpoint. The scientific archive contains checkpoints and optimizer histories, clean inputs, raw-source identities, committed predictions, tables, figure data and a resumable staged entry point. Executed stages, failed attempts and administrative costs are distinguished. No new wet-lab measurement was performed. The calibrated theoretical estimator is separately specified with complete proofs and arithmetic qualifications; no certified backend is claimed.
\section*{Ethics statement}
All measurements are from existing public experimental reports. No new participant, animal or laboratory experiment was performed. No clinical efficacy, biological mechanism, causal effect of an isolated chemical group, or theorem applicability is inferred from prediction errors. Public assay metadata and representation gaps are retained rather than filled with invented values.
\section*{Use of AI assistance}
AI assistance was used for code, analysis, manuscript drafting and evidence organization. Claims remain limited by the recorded source checks and executed computations. The final workflow illustration is reserved for separate preparation; this draft supplies its scientific specification.
''')
print(P)
