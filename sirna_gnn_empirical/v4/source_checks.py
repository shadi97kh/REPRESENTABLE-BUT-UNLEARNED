"""Bounded primary-source, citation and official-code update checks; no new baseline fit."""
from common import *
import requests,xml.etree.ElementTree as ET,datetime,re,configparser,urllib.parse,ast
V3=ROOT/'runs/sirna_gnn_empirical/v3-20260914T013314Z';out=RUN/'source_checks';out.mkdir(exist_ok=True)
queries={'iclr_author_guidelines':'https://iclr.cc/Conferences/2027/AuthorGuidelines','iclr_ai_policy':'https://iclr.cc/Conferences/2027/AIPolicyForAuthors','ensi_latest':'https://api.github.com/repos/tanwenchong/ENsiRNA/commits/HEAD','meg_latest':'https://api.github.com/repos/YuantingChen111/MEG-mod/commits/HEAD','davis_corrections':'https://api.crossref.org/works/10.1093/nar/gkaf479'}
access=[]
for name,url in queries.items():
 t=time.monotonic()
 try:
  r=requests.get(url,timeout=40);p=out/(name+('.json' if name.endswith('latest') or name=='davis_corrections' else '.html'));p.write_bytes(r.content);access.append(dict(name=name,url=url,status=r.status_code,bytes=len(r.content),sha256=sha(p),wall_s=time.monotonic()-t))
 except requests.RequestException as e:access.append(dict(name=name,url=url,error=str(e),wall_s=time.monotonic()-t))
write(out/'access.json',access)
# No known-failing command is rerun. Determine whether official code actually changed.
oldstatus=json.loads((V3/'sources/reproduction_status.json').read_text());changes=[]
for name,method in [('ensi_latest','ENsiRNA_mod'),('meg_latest','MEG_mod')]:
 p=out/(name+'.json');new=json.loads(p.read_text()).get('sha') if p.exists() else None;previous=oldstatus[method]['commit'];changes.append(dict(method=method,pinned_commit=previous,current_commit=new,official_change_detected=bool(new and new!=previous),previous_blocker=oldstatus[method]['blocker'],action='inspect changed official implementation before any fit' if new and new!=previous else 'no new official code supporting a resolved blocker; no repeated failed invocation'))
write(out/'published_update_checks.json',changes)
# Independently inspect all retained primary evidence bytes and current citation uses.
source=ROOT/'papers/interaction_recoverability_iclr2027/v9/source';tex=(source/'main.tex').read_text()+(source/'appendix.tex').read_text();bib=(source/'references.bib').read_text();keys=set(re.findall(r'@\w+\s*\{\s*([^,]+),',bib));used={k.strip() for z in re.findall(r'\\cite\w*(?:\[[^]]*\])?\{([^}]+)\}',tex) for k in z.split(',')};assert used<=keys
old=json.loads((V3/'citation_audit.json').read_text());audit=[]
for rec in old:
 for ev in rec['evidence']:assert sha(V3/ev['path'])==ev['sha256']
 texts=[]
 for ev in rec['evidence']:
  p=V3/ev['path']
  if p.suffix in ['.txt','.json']:texts.append(p.read_text(errors='replace'))
 text='\n'.join(texts);assert len(text)>50,rec['key']
 snippets=[re.sub(r'\s+',' ',m.group(0))[:350] for m in list(re.finditer(r'.{0,100}(?:siRNA|graph|normalization|regularization|regression|bootstrap|dropout|learning|prediction|neural|RNA|ensemble|optimizer|reproducib).{0,200}',text,re.I))[:3]]
 audit.append(dict(**rec,used_in_minimal_source=rec['key'] in used,retained_primary_bytes_verified=True,primary_context_excerpts=snippets,revision_scope='Same narrow contextual/software claims; primary Bramsen assay and source adjudication independently rechecked. No imported GNN/biological guarantee.'))
write(out/'citation_recheck.json',audit);assert keys=={x['key'] for x in audit}
# Retain the exact primary text establishing assay/replicate information and orientation limits.
for article,needle in [('PMC2685080',['triplicate','normalized','strand','10 nM']),('PMC12205987',['antisense','Supplementary Table S1','reverse complement'])]:
 root=ET.parse(V3/'sources'/f'{article}.xml').getroot();paras=[''.join(p.itertext()) for p in root.iter('p')];chosen=[p for p in paras if any(n.lower() in p.lower() for n in needle)];(out/f'{article}_methods_excerpts.txt').write_text('\n\n'.join(chosen))
# Retain provisional source limitations until exact released row/assay provenance resolves them.
write(out/'summary.json',dict(citation_entries=len(audit),cited_entries=len(used),new_published_fits=0,published_code_updates=changes,S1='Still quarantined: source text does not justify converting its target-sense printed antisense chemical strings.',primary_assay_replicates='Bramsen reports triplicate assays repeated twice; this does not recover endpoint covariance or independent-study replication.',iclr='Official 2027 pages checked live: nine-page main limit; references/appendix excluded; AI, ethics (max one page) and reproducibility statements explicitly excluded. Authors remain accountable for AI-assisted work.',scope='No private-anchor theorem or new wet-lab validation; no contact or publishing.'))
print(json.dumps(json.loads((out/'summary.json').read_text()),indent=2))
