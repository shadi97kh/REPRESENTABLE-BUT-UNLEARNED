import argparse,json,re,math,collections,itertools
from pathlib import Path
import pandas as pd
from common import *
ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);r=ap.parse_args().run
rows=[];excluded=[];rawrows=[];seen={};overlap={}
for filename in ['train_88_1.xlsx','valid_88_1.xlsx','test_88.xlsx','patent.xlsx']:
 p=r/'raw/ensi/ENsiRNA-mod/dataset'/filename;frame=pd.read_excel(p).fillna('');overlap[filename]=set()
 for i,series in frame.iterrows():
  x={str(k):v for k,v in series.to_dict().items()};rid=f'ensi:{filename}:{i+2}';rawrows.append({'record_id':rid,'source_row':x})
  overlap[filename].add((x['anti raw seq'],x['sense raw seq']))
  try:
   source=str(x.get('source','')).strip();assert re.fullmatch(r'\d{7,9}|US20120088815A1, EP2415869A1',source),'unresolved source identity'
   a=ensi_strand(x['anti raw seq'],x['anti mod'],x['anti pos']);s=ensi_strand(x['sense raw seq'],x['sense mod'],x['sense pos'])
   pairs,score=pairing(a['sequence'],s['sequence']);assert score>=.80,'unresolved duplex pairing'
   y=float(x['PCT']);dose=float(x['cc']);assert math.isfinite(y) and dose>0,'invalid outcome/concentration'
   obj={'dataset':'ENsiRNA','source_family':source,'record_id':rid,'source_id':str(x['ID']),'source_file':str(p.relative_to(r)),'source_row':int(i+2),'guide':a,'passenger':s,'pairing':pairs,'pairing_fraction':score,'outcome_raw':y,'activity':y/100,'sd_raw':None,'replicate_n':None,'dose_reported':dose,'dose_unit':'not_verified','target':None,'cell':None,'time_h':None,'delivery':None,'assay_id':None,'original_supplied_split':filename,'chemistry_resolved_at':'reported positional names; linkage direction and stereochemistry incomplete','source_table_verified':False}
   identity=digest([source,a,s,dose,y])
   if identity in seen:excluded.append({'record_id':rid,'reason':'exact duplicated molecular/dose/outcome record','duplicate_of':seen[identity]});continue
   seen[identity]=rid;rows.append(obj)
  except (AssertionError,ValueError,TypeError) as e:excluded.append({'record_id':rid,'reason':str(e)})
cm=pd.read_csv(r/'raw/cmsirnadb_APP.tsv',sep='\t',keep_default_na=False)
for i,series in cm.iterrows():
 x=series.to_dict();rid=f'cmAPP:{i+2}';rawrows.append({'record_id':rid,'source_row':x})
 try:
  assert x['Cell_Type']!='Transgenic mice','animal assay excluded from cell prediction'
  assert x['Time_of_administration']=='24h','unresolved incubation time'
  a=cm_strand(x['Antisense_seqence'],x['Modification_Types_Antisense_strand'],x['Modifications_AntiSense_strand']);s=cm_strand(x['Sense_seqence'],x['Modification_Types_Sense_strand'],x['Modifications_sense_strand'])
  pairs,score=pairing(a['sequence'],s['sequence']);assert score>=.8,'unresolved duplex pairing'
  y=float(x['Inhibition']);dose=float(x['Concentration'].removesuffix('nM'));assert math.isfinite(y) and dose>0,'invalid outcome/concentration'
  obj={'dataset':'CMsiRNAdb_APP','source_family':'Alnylam_APP_2020_2022_related','patent':x['patent_ID'],'record_id':rid,'source_id':x['ID'],'agent':x['The_name_of_double_helix'],'source_file':'raw/cmsirnadb_APP.tsv','source_row':int(i+2),'guide':a,'passenger':s,'pairing':pairs,'pairing_fraction':score,'outcome_raw':y,'activity':y/100,'sd_raw':float(x['SD']) if x['SD'] else None,'replicate_n':None,'dose_reported':dose,'dose_unit':'nM','target':x['Target_Gene'],'accession':x['Accession_number'],'cell':x['Cell_Type'],'time_h':24,'delivery':'RNAiMAX per source screening protocol; row/table linkage unresolved','assay_id':x['patent_ID']+'|'+x['Cell_Type']+'|'+x['Concentration']+'|24h','original_supplied_split':None,'chemistry_resolved_at':'indexed nucleotide descriptions and raw sequence agreement; unspecified stereochemistry retained','source_table_verified':False}
  identity=digest([obj['patent'],a,s,obj['cell'],dose,y,obj['agent']])
  if identity in seen:excluded.append({'record_id':rid,'reason':'exact duplicated molecular/assay/outcome record','duplicate_of':seen[identity]});continue
  seen[identity]=rid;rows.append(obj)
 except (AssertionError,ValueError,TypeError) as e:excluded.append({'record_id':rid,'reason':str(e)})
# Supporting Davis file is not silently repaired: all primary S1 reported antisense strings fail 13-mer complementarity.
d=pd.read_excel(r/'raw/davis2025_supplements/Supplemental Tables.xlsx',sheet_name='Supplemental Table 1',header=5).fillna('');dc=collections.Counter()
for i,x in d.iterrows():
 a=''.join(re.findall(r'\([mfr]([ACGU])\)',x.iloc[4]));s=''.join(re.findall(r'\([mfr]([ACGU])\)',x.iloc[5]));t=x.iloc[6];rc=''.join(COMP[c] for c in t[::-1]);dc['rows']+=1;dc['antisense_13mer_matches_target_sense']+=int(any(a[j:j+13] in t for j in range(len(a)-12)));dc['antisense_13mer_matches_target_antisense']+=int(any(a[j:j+13] in rc for j in range(len(a)-12)))
 excluded.append({'record_id':f'Davis:S1:{i+7}','reason':'reported strand/base orientation inconsistent with target and duplex; no unverified repair','compound':x.iloc[0]})
write_json(r/'davis_orientation_audit.json',dict(dc))
# Audit source-family and sequence links before choosing the split.
dsu=DSU(len(rows));source_first={};kmers={}
for i,x in enumerate(rows):
 if x['source_family'] in source_first:dsu.join(i,source_first[x['source_family']])
 else:source_first[x['source_family']]=i
 for st in ['guide','passenger']:
  seq=x[st]['sequence'].replace('T','U')
  for k in {seq[j:j+13] for j in range(len(seq)-12)}:
   if k in kmers:dsu.join(i,kmers[k])
   else:kmers[k]=i
components=collections.defaultdict(list)
for i,x in enumerate(rows):components[dsu.find(i)].append(i)
summary={'admitted_records':len(rows),'datasets':dict(collections.Counter(x['dataset'] for x in rows)),'sources':dict(collections.Counter(x['source_family'] for x in rows)),'exclusion_reasons':dict(collections.Counter(x['reason'] for x in excluded)),'supplied_ENsiRNA_raw_duplex_train_test_overlap':len(overlap['train_88_1.xlsx']&overlap['test_88.xlsx']),'global_study_sequence_components':[{'rows':len(ids),'source_families':sorted(set(rows[i]['source_family'] for i in ids)),'datasets':sorted(set(rows[i]['dataset'] for i in ids))} for ids in sorted(components.values(),key=len,reverse=True)],'B2_status':'primary experiment/table and chemistry comparability audit pending; provisional database blocks not admitted','B3_status':'no admitted four-condition factorial set','Davis_orientation':dict(dc)}
save_jsonl(r/'observations.jsonl',rows);save_jsonl(r/'raw_rows.jsonl',rawrows);save_jsonl(r/'eligibility_exclusions.jsonl',excluded);write_json(r/'data_identity.json',summary)
write_json(r/'reconcile.outputs.json',['observations.jsonl','raw_rows.jsonl','eligibility_exclusions.jsonl','data_identity.json','davis_orientation_audit.json'])
print(json.dumps(summary,indent=2))
