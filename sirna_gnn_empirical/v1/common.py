import hashlib,json,re
from pathlib import Path

def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
def write_json(p,x):p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,indent=2,allow_nan=False,ensure_ascii=False)+'\n')
def save_jsonl(p,x):Path(p).write_text(''.join(json.dumps(v,allow_nan=False,ensure_ascii=False)+'\n' for v in x))
def read_jsonl(p):return [json.loads(s) for s in Path(p).read_text().splitlines()]
COMP={'A':'U','C':'G','G':'C','U':'A','T':'A'}
def pairing(a,s):
 # Coordinates retain reported 5'-3' ordering. Pairing is a sequence-derived alignment, not a measured conformation.
 opts=[]
 for offset in range(-len(s)+1,len(a)):
  pairs=[(i,len(s)-1-(i-offset)) for i in range(len(a)) if 0<=i-offset<len(s)]
  if len(pairs)<14:continue
  matches=sum(COMP[a[i]]==s[j].replace('T','U') for i,j in pairs)
  opts.append((matches, -sum(COMP[a[i]]!=s[j].replace('T','U') for i,j in pairs),len(pairs),-abs(offset),-offset,pairs))
 if not opts:return [],0.
 best=max(opts);return best[-1],best[0]/best[2]
ALIASES={'2-O-Methyl ribose':'2-O-Methyl','2-Methoxy':'2-O-Methyl','2-Deoxy-2-Fluoro':'2-Fluoro','2-Deoxyribonucleotide':'2-Deoxy','2-Deoxythymidine':'2-Deoxy','UnLocked nucleic acid':'Unlocked nucleic acid','5-phosphate ribose':'5-Phosphate'}
def ensi_strand(seq,mods,positions):
 seq=str(seq).strip();assert re.fullmatch('[ACGUT]+',seq), 'noncanonical base string'
 n=len(seq);nodes=[{'base':b,'mods':[],'stereo':'not_reported'} for b in seq];link=['PO_assumed_from_annotation']*(n-1);terminal=[]
 if str(mods).strip() in ['0','0.0']:
  assert str(positions).strip() in ['0','0.0'],'inconsistent explicit-zero annotation'
 else:
  assert str(mods).strip() and str(positions).strip(),'missing chemical annotation'
  mm=[x.strip() for x in str(mods).split('*')];pp=[x.strip() for x in str(positions).split('*')];assert len(mm)==len(pp),'modification/position arity'
  for m,pos in zip(mm,pp):
   m=ALIASES.get(m,m)
   assert m not in ['Mutation','Inverted abasic','Assymetric siRNA (3-End Overhang)'],'unresolved base/geometry-changing annotation'
   for token in pos.split(','):
    i=int(float(token.strip()))-1;assert 0<=i<n,'modification index out of range'
    nodes[i]['mods'].append(m)
    if m in ['Phosphorothioate','Boranophosphate']:
     # Source does not resolve direction of indexed linkage: retain positional tag, do not invent a resolved bond.
     nodes[i]['linkage_direction_unresolved']=True
    if m in ['5-Phosphate','Dodecyl derivative']:terminal.append([i,m])
 for node in nodes:
  node['mods']=sorted(set(node['mods'])) or ['unmodified_as_annotated']
 return {'sequence':seq,'nodes':nodes,'linkages':link,'terminal':terminal,'annotation_scope':'Positional modification names from released workbook; unstated chemistry/stereochemistry is not resolved.'}
BASE_WORDS={'adenosine':'A','guanosine':'G','cytidine':'C','uridine':'U','thymidine':'T'}
def cm_strand(seq,description,raw):
 seq=str(seq).strip();assert re.fullmatch('[ACGUT]+',seq),'noncanonical base string';assert raw and description,'missing chemical annotation'
 nodes=[None]*len(seq);link=['PO']*(len(seq)-1);terminal=[]
 for token in description.split(' || '):
  if 'GalNAc' in token and not re.match(r'^\d+\*',token):
   assert raw.endswith('L96'),'unlocalized conjugate'
   terminal.append([len(seq)-1,'L96_GalNAc']);continue
  m=re.fullmatch(r'(\d+)\*(.+)',token.strip());assert m,'unparsed indexed chemistry description'
  idx=int(m[1])-1;desc=m[2]
  if 'GalNAc' in desc:
   terminal.append([idx,'L96_GalNAc']);continue
  base=next((v for k,v in BASE_WORDS.items() if k in desc),None);assert base and 0<=idx<len(seq),'unknown chemistry/base index';assert seq[idx]==base,'description/base disagreement'
  if 'hexadecyl' in desc:mod='2-O-hexadecyl'
  elif 'Glycol Nucleic Acid' in desc:mod='Glycol nucleic acid'
  elif 'Deoxy' in desc:mod='2-Deoxy'
  elif 'Methyl' in desc:mod='2-O-Methyl'
  elif 'Fluoro' in desc:mod='2-Fluoro'
  else:raise ValueError('unresolved sugar chemistry: '+desc)
  ms=[mod]
  if 'Phosphorothioate' in desc:
   ms.append('Phosphorothioate')
   if idx<len(link):link[idx]='PS'
   else:terminal.append([idx,'terminal_PS'])
  if 'Vinyl' in desc:ms.append('5-vinylphosphonate');terminal.append([idx,'5-vinylphosphonate'])
  assert nodes[idx] is None,'duplicate nucleotide index'
  nodes[idx]={'base':base,'mods':ms,'stereo':'S_GNA' if base=='T' and mod=='Glycol nucleic acid' and '(Tgn)' in raw else 'not_reported'}
 assert all(x is not None for x in nodes),'missing nucleotide chemical description'
 return {'sequence':seq,'nodes':nodes,'linkages':link,'terminal':terminal,'annotation_scope':'Database indexed description checked against bases; original case-sensitive notation preserved. Unspecified stereochemistry remains unresolved.'}
class DSU:
 def __init__(self,n):self.p=list(range(n))
 def find(self,x):
  while self.p[x]!=x:self.p[x]=self.p[self.p[x]];x=self.p[x]
  return x
 def join(self,a,b):self.p[self.find(a)]=self.find(b)
