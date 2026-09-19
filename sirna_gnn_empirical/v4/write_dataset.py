from common import *
import copy
out=RUN/'dataset';out.mkdir(exist_ok=True);targets=pd.read_csv(RUN/'adjudication/primary_reanchored_targets.csv').set_index('record_id').to_dict('index');admitted=[];excluded=[]
for i,x in enumerate(rows_all()):
 if x['dataset']!='ENsiRNA':continue
 if x['source_family']=='19282453' and x['record_id'] not in targets:excluded.append(dict(record_id=x['record_id'],reason='No uniquely resolved primary measurement for reported chemical identity'));continue
 z=copy.deepcopy(x);z['original_activity']=z['activity'];z['target_version']='primary_assay_sensitivity_v1'
 if x['record_id'] in targets:z['activity']=float(targets[x['record_id']]['primary_activity']);z['target_provenance']=targets[x['record_id']]
 else:z['target_provenance']='unchanged released label; not newly adjudicated'
 z['original_graph_index']=i;admitted.append(z)
assert len(admitted)==2607 and len(excluded)==320
p=out/'primary_reanchored_observations.jsonl';p.write_text(''.join(json.dumps(x)+'\n' for x in admitted));write(out/'quarantined_source_rows.json',excluded);write(out/'manifest.json',dict(version='primary_assay_sensitivity_v1',rows=len(admitted),Bramsen_rows=len(targets),other_source_rows=1003,excluded_Bramsen=320,sha256=sha(p),target_map_sha256=sha(RUN/'adjudication/primary_reanchored_targets.csv'),scope='A primary-assay sensitivity dataset with retained original values; not a certified correction to the released benchmark.',original_data_unchanged=True));print('Versioned dataset:',len(admitted),'rows; original values retained.')
