"""Primary-only Bramsen rectangle identity reconciliation; no fitting or predicted labels."""
from common import *
import xml.etree.ElementTree as E,collections
out=RUN/'b3';out.mkdir(exist_ok=True);source=RUN/'sources/PMC2685080';root=E.fromstring((RUN/'sources/PMC2685080.xml').read_text());table=root.find('.//table-wrap[@id="T1"]')
# Exact named chemical aliases from primary Table1 footnote to the frozen positional vocabulary.
alias={'OMe':'2-O-Methyl','F':'2-Fluoro','DNA':'2-Deoxy','AEM':'2-Aminoethoxymethyl','APM':'2-Aminopropoxymethyl','EA':'2-Aminoethyl','CE':'2-Cyanoethyl','GE':'2-Guanidinoethyl','HM':'4-C-Hydroxymethyl Deoxyribonucleic acid','LNA':'Locked nucleic acid','ALN':'Alfa-L-Locked nucleic acid','ADA':'2-N-Adamant-1-yl-Methylcarbonyl-2-Amino-Locked nucleic acid','PYR':'2-N-Pyren-1-yl Methyl-2-Amino-Locked nucleic acid','OX':'Oxetane-Locked nucleic acid','CENA':'2,4-Carbocyclic-Ethylene-bridged nucleic acid-Locked nucleic acid','CLNA':'2,4-Carbocyclic-Locked nucleic acid-Locked nucleic acid','UNA':'Unlocked nucleic acid','ANA':'Altritol nucleic acid','AENA':'2-Deoxy-2-N,4-C-Ethylene-Locked nucleic acid','HNA':'Hexitol nucleic acid'}
# AP has sugar/base ambiguity in the retained source vocabulary and is deliberately not guessed.
write(out/'mapping_rules.json',dict(aliases=alias,AP='unresolved sugar/base mapping excluded',rounding_tolerance_activity=.000500001,reference_AS='W053',reference_SS='W207',exact_bases=True,source_scale='Relative eGFP levels are fractions, not percentages; activity=1-value',selection='Identity first; rounding check validates primary/database reconciliation, never prediction fit quality'))
strands={};ex=[]
for tr in table.findall('.//tr'):
 td=tr.findall('td')
 if len(td)!=6:continue
 name=''.join(td[0].itertext()).strip();strand=''.join(td[1].itertext()).strip();nodes=[];unresolved=[]
 def addchars(s):
  for c in (s or '').replace(' ','').replace('\n',''):
   if c in 'ACGUT':nodes.append(dict(base=c,mods=[]))
   elif c not in ['′',"'",'5','3','–','-']:unresolved.append(c)
 def walk(el):
  addchars(el.text)
  for child in el:
   if child.tag in ['sub','sup']:
    chem=''.join(child.itertext()).strip()
    if not nodes or chem not in alias:unresolved.append('chemical:'+chem)
    else:nodes[-1]['mods'].append(alias[chem])
   else:walk(child)
   addchars(child.tail)
 walk(td[5]);seq=''.join(n['base'] for n in nodes)
 if unresolved:ex.append(dict(name=name,reason='primary markup/chemical alias unresolved',details=','.join(unresolved)));continue
 # A bold nucleotide without an explicit chemical subscript is not silently treated as unmodified.
 # Current T1 entries tag modifications using subscript; all parsed states are retained with XML evidence.
 strands[name]=dict(name=name,strand=strand,sequence=seq,nodes=nodes,xml=E.tostring(td[5],encoding='unicode'))
write(out/'primary_strand_states.json',strands)
def key(st):return (st['sequence'],tuple(tuple(sorted(m for m in n['mods'] if m!='unmodified_as_annotated')) for n in st['nodes']))
for name in ['W053','W207']:assert name in strands
am={};sm={}
for name,st in strands.items():
 target=am if st['strand']=='AS' else sm;target.setdefault(key(st),[]).append(name)
rows=[x for x in rows_all() if x['source_family']=='19282453'];mapping=collections.defaultdict(list);rowmap=[]
for x in rows:
 aa=am.get(key(x['guide']),[]);ss=sm.get(key(x['passenger']),[])
 if len(aa)==1 and len(ss)==1:mapping[(aa[0],ss[0])].append(x)
 rowmap.append(dict(record_id=x['record_id'],antisense_candidates='|'.join(aa),sense_candidates='|'.join(ss),identity_unique=len(aa)==1 and len(ss)==1))
export(pd.DataFrame(rowmap),'b3/record_identity_mapping.csv')
df=pd.read_excel(source/'gkp106_nar-02440-c-2008-File010.xls').iloc[:,:7];df.columns=['AS','SS','pair','remaining','sd','viability','viability_sd'];primary={};reconciled={};audit=[]
for (a,s),g in df.groupby(['AS','SS']):
 if len(g)!=1:ex.append(dict(name=a+'|'+s,reason='duplicate primary pair identity',details=len(g)));continue
 r=g.iloc[0];obs=mapping.get((a,s),[])
 if len(obs)!=1:
  audit.append(dict(AS=a,SS=s,admitted_records=len(obs),status='ambiguous_or_missing_admitted_identity'));continue
 x=obs[0];y=1-float(r.remaining);err=abs(y-x['activity']);ok=np.isfinite(y) and err<=.000500001
 audit.append(dict(AS=a,SS=s,record_id=x['record_id'],primary_activity=y,archived_activity=x['activity'],absolute_difference=err,status='verified' if ok else 'source_mean_mismatch'))
 if ok:reconciled[(a,s)]=dict(record_id=x['record_id'],activity=y,measurement_sd=float(r.sd),AS=a,SS=s,sequence_group=x['sequence_group'],study_group=x['study_group'])
export(pd.DataFrame(audit),'b3/primary_reconciliation.csv');export(pd.DataFrame(ex),'b3/exclusions.csv')
ref=('W053','W207');rect=[]
for a in sorted(strands):
 for s in sorted(strands):
  if a=='W053' or s=='W207' or strands[a]['strand']!='AS' or strands[s]['strand']!='SS':continue
  if strands[a]['sequence']!=strands['W053']['sequence'] or strands[s]['sequence']!=strands['W207']['sequence']:continue
  keys=[ref,(a,'W207'),('W053',s),(a,s)]
  if not all(k in reconciled for k in keys):continue
  v=[reconciled[k] for k in keys];assert len({z['record_id'] for z in v})==4
  rect.append(dict(rectangle_id=hashlib.sha256((a+'|'+s).encode()).hexdigest()[:16],AS=a,SS=s,record_ids={k:z['record_id'] for k,z in zip(['00','10','01','11'],v)},observed={k:z['activity'] for k,z in zip(['00','10','01','11'],v)},endpoint_sd={k:z['measurement_sd'] for k,z in zip(['00','10','01','11'],v)},observed_interaction=v[3]['activity']-v[1]['activity']-v[2]['activity']+v[0]['activity'],sequence_group=v[0]['sequence_group'],study_group=v[0]['study_group']))
write(out/'admitted_rectangles.json',rect);write(out/'admission_summary.json',dict(primary_rows=len(df),parsed_primary_strands=len(strands),archived_source_rows=len(rows),unique_mapped_records=sum(x['identity_unique'] for x in rowmap),reconciled_primary_pairs=len(reconciled),reference_reconciled=ref in reconciled,admitted_rectangles=len(rect),exact_sequence_backgrounds=len({z['sequence_group'] for z in rect}),scope='Retrospective primary-reconciled one-source/one-sequence panel; no independent rectangle intervals or causal/new-wetlab claim; original APP B3 unchanged',protocol_hash=sha(RUN/'b3_protocol_addendum.json')))
print(json.loads((out/'admission_summary.json').read_text()))
