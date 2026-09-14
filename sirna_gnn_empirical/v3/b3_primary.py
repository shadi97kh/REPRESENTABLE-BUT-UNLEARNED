"""Standalone primary measured panel; archived release remains unchanged."""
from common import *
import collections,datetime
out=RUN/'b3_primary';out.mkdir(exist_ok=True)
# Explicit source-driven amendment: prior admission required an archived reference that is absent.
protocol=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),supersedes_for_this_new_primary_panel='b3_protocol_addendum.json restricted to already admitted four endpoints; that attempt remains archived and has zero admitted rectangles.',reason='Primary W053/W207 reference and complete measured combinations exist but the unmodified reference is absent from the chemistry-only benchmark release. Requiring that preexisting row is not a biological admission requirement. Admit independently verified primary states rather than fabricate an archived identity.',source='Bramsen2009 primary Table1 states and SupplementalTable1 measurements',source_hash=sha(RUN/'sources/PMC2685080/gkp106_nar-02440-c-2008-File010.xls'),state_hash=sha(RUN/'b3/primary_strand_states.json'),selection='Every parseable AS and SS with exact base sequence/length of W053/W207, no unclear chemical aliases, and all four finite means. No measured effect magnitude or prediction error selects a rectangle.',label_scale='1 - primary relative eGFP level; source SD unchanged in fractional units; no clipping',known_discrepancy='Primary/released outcome mismatch is reported separately; no existing campaign label is overwritten.',representation='Source-established positional names and bases in frozen 93-column vocabulary; annotation-only unresolved backbone, as original ENsiRNA representation, because primary table does not specify bond directions. Unmodified node names explicitly set; original deterministic duplex pairing.',models='Only newly frozen outer0 models whose training/validation exclude source19282453 and every primary-panel sequence 13mer; no refit or tuning. Use checkpoint preprocessing unchanged.',metrics=['interaction MSE/MAE/correlation/spread/bias','zero interaction reference','ten seed risks and ensemble risk','all endpoint and interaction predictions'],uncertainty='One source and one sequence background with shared controls. No independent-rectangle bootstrap, no per-interaction SE without covariance.',scope='Retrospective measured four-condition assay-scale contrasts; not isolated GNA, biological mechanism, private-anchor theorem applicability or new wet-lab validation.')
if (out/'protocol.json').exists():protocol=json.loads((out/'protocol.json').read_text())
else:write(out/'protocol.json',protocol)
state=json.loads((RUN/'b3/primary_strand_states.json').read_text());df=pd.read_excel(RUN/'sources/PMC2685080/gkp106_nar-02440-c-2008-File010.xls').iloc[:,:7];df.columns=['AS','SS','pair','remaining','sd','viability','viability_sd'];values={}
for key,g in df.groupby(['AS','SS']):
 if len(g)==1 and np.isfinite(g.iloc[0].remaining):values[key]=g.iloc[0]
aset=[a for a,x in state.items() if x['strand']=='AS' and x['sequence']==state['W053']['sequence']];sset=[s for s,x in state.items() if x['strand']=='SS' and x['sequence']==state['W207']['sequence']];rect=[];endpoint=set()
for a in sorted(aset):
 for s in sorted(sset):
  if a=='W053' or s=='W207':continue
  keys=[('W053','W207'),(a,'W207'),('W053',s),(a,s)]
  if not all(k in values for k in keys):continue
  # Only chemically distinct interventions, not duplicated names for an identical state.
  if state[a]['nodes']==state['W053']['nodes'] or state[s]['nodes']==state['W207']['nodes']:continue
  val=[1-float(values[k].remaining) for k in keys];ids=['bramsen:'+k[0]+'|'+k[1] for k in keys];endpoint.update(keys);rect.append(dict(rectangle_id=hashlib.sha256((a+'|'+s).encode()).hexdigest()[:16],AS=a,SS=s,record_ids=dict(zip(['00','10','01','11'],ids)),observed=dict(zip(['00','10','01','11'],val)),endpoint_sd={lab:float(values[k].sd) for lab,k in zip(['00','10','01','11'],keys)},observed_interaction=val[3]-val[1]-val[2]+val[0]))
# Graph generation is label-independent and uses only exact primary states and the fixed schema.
manifest=json.loads((OLD/'graph_manifest.json').read_text());vocab=manifest['modification_vocabulary'];X=np.zeros((len(endpoint),64,93),np.float32);A=np.zeros((len(endpoint),8,64,64),np.float32);mask=np.zeros((len(endpoint),64),np.float32);observations=[];comp={'A':'U','C':'G','G':'C','U':'A','T':'A'}
def align(g,s):
 candidates=[]
 for off in range(-len(s)+1,len(g)):
  pp=[(i,len(s)-1-(i-off)) for i in range(len(g)) if 0<=i-off<len(s)]
  if len(pp)<14:continue
  matches=sum(comp[g[i]]==s[j].replace('T','U') for i,j in pp);candidates.append(((matches,-(len(pp)-matches),len(pp),-abs(off),-off),pp))
 rank,pp=max(candidates);assert rank[0]/len(pp)>=.8;return pp,rank[0]/len(pp)
for b,(a,s) in enumerate(sorted(endpoint)):
 row=dict(record_id='bramsen:'+a+'|'+s,dataset='Bramsen_primary',source_family='19282453',study_group='Bramsen_primary',sequence_group='Bramsen_primary_single_background',activity=1-float(values[a,s].remaining),measurement_sd=float(values[a,s].sd),dose_reported=10.,dose_unit='nM',time_h=72,cell='HeLa-eGFP',assay_id='Bramsen2009_eGFP_10nM_72h',delivery='INTERFERin',AS=a,SS=s)
 for st,key,n in [(0,'guide',a),(1,'passenger',s)]:
  item=state[n];nodes=[]
  for i,node in enumerate(item['nodes']):
   mods=node['mods'] or ['unmodified_as_annotated'];assert all(m in vocab for m in mods);nodes.append(dict(base=node['base'],mods=mods,stereo='not_reported'));j=32*st+i;X[b,j,'ACGUT'.index(node['base'])]=1;X[b,j,5+st]=1;X[b,j,7+i]=1;mask[b,j]=1
   for m in mods:X[b,j,39+vocab.index(m)]=1
   X[b,j,84]=1;X[b,j,87]=i==0;X[b,j,88]=i==len(item['nodes'])-1
   if i<len(item['nodes'])-1:A[b,2,j+1,j]=1;A[b,5,j,j+1]=1
  row[key]=dict(sequence=item['sequence'],nodes=nodes,linkages=['PO_assumed_from_annotation']*(len(nodes)-1),terminal=[])
 pp,fr=align(row['guide']['sequence'],row['passenger']['sequence']);row['pairing']=pp;row['pairing_fraction']=fr
 for i,j in pp:
  rel=6 if comp[row['guide']['sequence'][i]]==row['passenger']['sequence'][j].replace('T','U') else 7;A[b,rel,i,32+j]=A[b,rel,32+j,i]=1
 observations.append(row)
(out/'observations.jsonl').write_text(''.join(json.dumps(x)+'\n' for x in observations));np.savez_compressed(out/'graphs.npz',X=X,A=A,mask=mask,C=np.zeros((len(observations),8),np.float32));write(out/'rectangles.json',rect)
write(out/'admission_summary.json',dict(measured_conditions=len(observations),complete_rectangles=len(rect),eligible_AS=len(aset),eligible_SS=len(sset),exact_sequence_backgrounds=1,source_families=1,protocol_hash=sha(out/'protocol.json'),original_APP_rectangles=0,proof_of_measurement='Every endpoint comes directly from the primary matrix; no model output is used as a label.'))
print(json.loads((out/'admission_summary.json').read_text()))
