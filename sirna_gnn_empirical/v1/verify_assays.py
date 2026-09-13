"""Reconcile APP numeric outcomes and chemical strings against the related primary tables."""
import argparse,json,re,collections,itertools,math
from pathlib import Path
from common import *
ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);r=ap.parse_args().run
obs=read_jsonl(r/'observations.jsonl');raw={x['record_id']:x['source_row'] for x in read_jsonl(r/'raw_rows.jsonl')}
p20=json.loads((r/'source_inspection/US20220002724A1_tables.json').read_text());p22=json.loads((r/'source_inspection/US20240117349A1_tables.json').read_text())
lookup={}
def add(agent,cell,dose,mean,sd,pub,table,row):
 try:lookup[(pub,agent,cell,float(dose))]={'mean_remaining':float(mean),'sd_remaining':float(sd),'primary_publication':pub,'table':table,'table_row':row}
 except ValueError:pass
for ti,num,cells in [(53,'4',['Primary Cynomolgus Monkey Hepatocytes','Be(2)C cell line']),(64,'7',['Primary mouse hepatocytes','Neuro2A cell line'])]:
 for rownum,a in enumerate(p20[ti]['rows']):
  if len(a)==9 and a[0].startswith('AD-'):
   for ci,cell in enumerate(cells):
    for di,dose in enumerate([10,.1]):
     k=1+4*ci+2*di;add(a[0],cell,dose,a[k],a[k+1],'US20220002724A1',num,rownum)
for ti,num,cell in [(90,'17','Be(2)C cell line'),(93,'18','Neuro2A cell line')]:
 for rownum,a in enumerate(p20[ti]['rows']):
  if len(a)==5 and a[0].startswith('AD-'):add(a[0],cell,a[3],a[1],a[2],'US20220002724A1',num,rownum)
for rownum,a in enumerate(p20[120]['rows']):
 if len(a)==5 and a[0].startswith('AD-'):
  for di,dose in enumerate([10,.1]):add(a[0],'Be(2)C cell line',dose,a[1+2*di],a[2+2*di],'US20220002724A1','29',rownum)
for rownum,a in enumerate(p22[15]['rows']):
 if len(a)==8 and a[0].startswith('AD-'):
  for di,dose in enumerate([10,1,.1]):add(a[0],'Be(2)C cell line',dose,a[2+2*di],a[3+2*di],'US20240117349A1','6',rownum)
def norm(s):return re.sub(r'\s+','',s)
def logical_rows(rows):
 current=[]
 for row in rows+[['']]:
  if not any(row):
   if current:yield [''.join(v) for v in zip(*current)];current=[]
  elif current and len(row)!=len(current[0]):
   yield [''.join(v) for v in zip(*current)];current=[row]
  else:current.append(row)
chem={}
for ti,ai,si,gi in [(44,0,1,3),(55,0,1,3),(73,0,1,3),(77,0,1,3),(81,1,2,4)]:
 for row in logical_rows(p20[ti]['rows']):
  if len(row)>max(ai,si,gi) and row[ai].startswith('AD-'):chem[('US20220002724A1',norm(row[ai]))]=(norm(row[gi]),norm(row[si]),ti)
for ti in [5,9]:
 for row in logical_rows(p22[ti]['rows']):
  if len(row)>=5 and row[0].startswith('AD-'):chem[('US20240117349A1',norm(row[0]))]=(norm(row[3]),norm(row[1]),ti)
accepted=[];checks=[];excluded=[];corrections=[]
for x in obs:
 if x['dataset']=='ENsiRNA':accepted.append(x);continue
 pub='US20220002724A1' if x['patent']=='WO2020132227A2' else 'US20240117349A1';key=(pub,x['agent'],x['cell'],x['dose_reported']);v=lookup.get(key)
 if not v:excluded.append({'record_id':x['record_id'],'reason':'no exact agent/cell/dose primary-table match'});continue
 primary_chem=chem.get((pub,x['agent']));a=raw[x['record_id']]
 match=bool(primary_chem and primary_chem[:2]==(norm(a['Modifications_AntiSense_strand']),norm(a['Modifications_sense_strand'])))
 check={'record_id':x['record_id'],**v,'chemical_string_match':match,'database_inhibition':x['outcome_raw'],'primary_inhibition':100-v['mean_remaining'],'database_sd':x['sd_raw']}
 checks.append(check)
 if not match:excluded.append({'record_id':x['record_id'],'reason':'primary chemical strings not exactly reconciled','primary_table':v['table']});continue
 newy=100-v['mean_remaining'];delta=newy-x['outcome_raw']
 if abs(delta)>.015:corrections.append({**check,'difference':delta,'rule':'Primary table supersedes database transformed outcome for an exactly matched molecule and assay.'})
 x.update(database_outcome_raw=x['outcome_raw'],outcome_raw=newy,activity=newy/100,sd_raw=v['sd_remaining'],source_table_verified=True,primary_evidence=v,assay_id=pub+'|Table'+v['table']+'|'+x['cell']+'|'+str(x['dose_reported'])+'nM|24h',delivery='RNAiMAX; same named primary screening protocol',replicate_n=None)
 x['experiment_id']=pub+'|Table'+v['table']+'|'+x['cell'];x['response_definition']='100 minus percent mRNA remaining relative to AD-1955 non-targeting control (2020 family), or source stated normalized mRNA (2022 family)';accepted.append(x)
# Loss-independent chemical contrasts: exact raw bases, same table/cell/dose/time/delivery and conjugates.
blocks=collections.defaultdict(list)
for x in accepted:
 if x['dataset']=='CMsiRNAdb_APP':
  k=(x['assay_id'],x['guide']['sequence'],x['passenger']['sequence'],digest([x['guide']['terminal'],x['passenger']['terminal']]))
  blocks[k].append(x)
pairs=[];rectangles=[]
def cv(x):
 return tuple((strand,i,tuple(n['mods']),n['stereo']) for strand in ['guide','passenger'] for i,n in enumerate(x[strand]['nodes']))
for k,group in sorted(blocks.items()):
 unique={cv(x):x for x in sorted(group,key=lambda x:x['record_id'])}
 vals=list(unique.items())
 for (va,a),(vb,b) in itertools.combinations(vals,2):
  changed=[{'strand':i[0],'position_1based':i[1]+1,'from':i[2],'to':j[2],'from_stereo':i[3],'to_stereo':j[3]} for i,j in zip(va,vb) if i!=j]
  if not changed:continue
  pairs.append({'pair_id':digest([a['record_id'],b['record_id']])[:16],'reference_id':a['record_id'],'changed_id':b['record_id'],'assay_id':a['assay_id'],'background_id':digest(k[1:])[:16],'changes':changed,'observed_difference':b['activity']-a['activity'],'measurement_sd_reference':a['sd_raw']/100,'measurement_sd_changed':b['sd_raw']/100,'mean_difference_se':None,'uncertainty_gap':'Replicate number and shared-control covariance not resolved; SDs are not SEs.'})
 # A and B each change one distinct nucleotide chemistry state; the exact union must be observed.
 for v00,x00 in vals:
  singles=[(v,x,[i for i,(a,b) in enumerate(zip(v00,v)) if a!=b]) for v,x in vals if v!=v00]
  singles=[z for z in singles if len(z[2])==1]
  for (va,xa,ia),(vb,xb,ib) in itertools.combinations(singles,2):
   if ia==ib:continue
   v11=list(v00);v11[ia[0]]=va[ia[0]];v11[ib[0]]=vb[ib[0]];x11=unique.get(tuple(v11))
   if x11:rectangles.append({'ids':[x00['record_id'],xa['record_id'],xb['record_id'],x11['record_id']],'assay_id':x00['assay_id'],'background_id':digest(k[1:])[:16]})
rectangles=list({tuple(sorted(x['ids'])):x for x in rectangles}.values())
save_jsonl(r/'verified_observations.jsonl',accepted);save_jsonl(r/'primary_row_checks.jsonl',checks);save_jsonl(r/'primary_exclusions.jsonl',excluded);save_jsonl(r/'primary_label_corrections.jsonl',corrections);save_jsonl(r/'eligible_pairs.jsonl',pairs);save_jsonl(r/'eligible_rectangles.jsonl',rectangles)
summary={'admitted':len(accepted),'datasets':dict(collections.Counter(x['dataset'] for x in accepted)),'verified_APP_assay_pools':dict(collections.Counter(x['assay_id'] for x in accepted if x['dataset']=='CMsiRNAdb_APP')),'exclusions':dict(collections.Counter(x['reason'] for x in excluded)),'primary_label_corrections':len(corrections),'chemical_pairs':len(pairs),'pair_backgrounds':len(set(x['background_id'] for x in pairs)),'two_single_position_rectangles':len(rectangles),'rectangle_backgrounds':len(set(x['background_id'] for x in rectangles))}
write_json(r/'verified_eligibility.json',summary);write_json(r/'verify_assays.outputs.json',['verified_observations.jsonl','primary_row_checks.jsonl','primary_exclusions.jsonl','primary_label_corrections.jsonl','eligible_pairs.jsonl','eligible_rectangles.jsonl','verified_eligibility.json']);print(json.dumps(summary,indent=2))
